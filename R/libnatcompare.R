##########################
#
# Functions (put these in a package)
#
##########################
getCbrfcWebDataUrl <- function(site.names, start = "1900-01-01", end = as.character(Sys.Date()))
{
		
	url <- paste("http://www.nwrfc.noaa.gov/westernwater/inc/getcsv.php?id=", site.names,
		"&tablelist=obs&dataobs=", start, "&dataobsend=", end, 
		"&datacharcs=&datacharcsend=&datafcst=&datafcstend=&dataens_memb=&dataens_membend=&dataobs_mean=&dataobs_meanend=&datacli_index=&datacli_indexend=&dataobsind=&ver=table&page=data", 
		sep='')
	return(url)
		
}

getCbrfcForecastDataUrl <- function(site.names, start = "1900-01-01", end = as.character(Sys.Date()))
{

n <- length(site.names)
script <- "http://www.nwrfc.noaa.gov/westernwater/inc/getcsv.php"

keys <- c('id', 'tablelist', 'dataobs', 'dataobsend', 'datacharcs',
'datacharcsend', 'datafcst', 'datafcstend', 'dataens_memb', 'dataens_membend', 
'dataobs_mean', 'dataobs_meanend', 'datacli_index',
'datacli_indexend', 'dataobsind', 'ver', 'page')

vals <- c('fcst', start, end, '1971', '', '1979-04-01', '2009-07-01', 
'2007-02-20', '2009-07-28', 'January', 'December', '1871-01-01', '2008-12-01', 
'MEI', '', 'data')

x <- character()
for(site in site.names){
	args <- paste(keys,c(site,vals),sep='=',collapse='&')
	x <- c(x,paste(script,args,sep='?'))
}

return(x)
#id=GBRW4&tablelist=GBRW4&dataobs=GBRW4&dataobsend=GBRW4&datacharcs=GBRW4&datacharcsend=GBRW4&datafcst=fcst&datafcstend=1900-01-01&dataens_memb=2009-10-07&dataens_membend=1971&dataobs_mean=&dataobs_meanend=1979-04-01&datacli_index=2009-07-01&datacli_indexend=2007-02-20&dataobsind=2009-07-28&ver=January&page=December&id=1871-01-01&tablelist=2008-12-01&dataobs=MEI&dataobsend=&datacharcs=data
#id=GBRW4&tablelist=fcst&dataobs=1957-01-01&dataobsend=2009-12-01&datacharcs=1971&datacharcsend=&datafcst=1979-04-01&datafcstend=2009-07-01&dataens_memb=2007-02-20&dataens_membend=2009-07-28&dataobs_mean=January&dataobs_meanend=December&datacli_index=1871-01-01&datacli_indexend=2008-12-01&dataobsind=MEI&ver=&page=data
}



downloadAndReadCbrfcWebData <- function(site.names, 
	download.dir = 'downloads', from.cache = TRUE, missing.val = -9999, ...)
{
	# This function takes a a character vector of river forcast center site 
	# names,and downloads the unregulated flow data for that site from the 
	# RFC. The data is read in and returned in a list
	
	urls <- getCbrfcWebDataUrl( site.names[!is.na(site.names)], ...)
	
	#Create dir of downloaded data
	if(!file.exists(download.dir)) dir.create( download.dir )

	#Create list of files
	files <- character()
	for(i in 1:length(site.names)){
		if(!is.na(site.names[i]))
			files <- c( files, file.path( download.dir, paste( site.names[i], 'csv',sep='.' )))
	}
	
	#Download and read in files
	data.cbrfc <- list()
	for(i in 1:length(urls)) {
		if( !from.cache || !file.exists(files[i]) )
			download.file( urls[i], destfile=files[i], 
				method = "wget", quiet = T)
		data.cbrfc[[i]] <- unique( read.csv( files[i] ) )
	}
	
	names(data.cbrfc) <- site.names[!is.na(site.names)]
	
	for(i in 1:length(data.cbrfc)){
			#data is listed with most recent first so get the last entry
		n <- length(as.character(data.cbrfc[[i]]$obsdate))
			#Get the start year and month
		s <- strsplit(as.character(data.cbrfc[[i]]$obsdate[n]),'-')[[1]]
		s <- as.integer(s)
			#must reverse the data for the same reason
		data.cbrfc[[i]] <- ts(rev(data.cbrfc[[i]]$val), 
								start=c(s[1], s[2]), frequency=12)
		data.cbrfc[[i]][ data.cbrfc[[i]] == missing.val ] = NA
	}
	
	# coerce time series list into multi dimensional time series
	
	# first find the time series with the largest range 
	widest <- 1
	r.widest <- range(time(data.cbrfc[[1]]))
	r.widest <- r.widest[2] - r.widest[1]
	for(i in 2:length(data.cbrfc)){
		r <- range(time(data.cbrfc[[i]]))
		r <- r[2] - r[1]
		if(r > r.widest){
			widest <- i
			s.widest <- start(data.cbrfc[[widest]])
		}
	
	}
	
	#initilize twice so range is range of longest ts
	mts <- cbind(data.cbrfc[[1]],data.cbrfc[[widest]])
	mts <- cbind(data.cbrfc[[1]],data.cbrfc[[2]])
	for(i in 2:(length(data.cbrfc)-1))		
		mts <- cbind(mts,data.cbrfc[[i+1]])
		
	mts <- ts(mts,start = s.widest, frequency = 12)
	attr(mts,'dimnames') <- list(NULL,site.names)	
	
	return(mts)
	
}

cumsum.ts <- function(x){
	
	# does cumsum on a time series, converts NA to zero
	# because NA will corrupt the remaining values in the sum 
	
	x[is.na(x)] = 0
	if(is.null(dim(x)))
			# Single dimension time series
		x <- ts(cumsum(as.vector(x)),start=start(x),frequency=frequency(x))
	else
			# Multi dimension time series
		for(i in 1:dim(x)[2])
			x[,i] <- cumsum(as.vector(x[,i]))
	return(x)
	
	
}

diff.ts <- function(x, y, cumsum = FALSE, percent = FALSE){
	
	# takes either single or multidimensional ts and performs 
	# differences, can return the strait difference, cumulative 
	# difference or percent difference 
	if( cumsum )
		z <- cumsum.ts(x-y)
	else if( percent ) 
		z <- 200*abs(x - y)/(x + y)
	else 
		z <- x - y
		
	attr(z,"dimnames") <- attributes(x)$dimnames
	return(z)
	
}

	# calculate a statistic for each period in a time series
	# default period id frequency but any window will work
seasonal.stat <- function(x, fun, win = frequency(x), ...){
	
	#figure out if ts is single or multi dimensional
	if(is.null(dim(x))){
			#Single ts
		return(seasonal.stat.work(x,fun=fun,win=win,giveStack=F,...))
	}else{
		
			# Multidimensional ts
			# get the seasonal stat for each ts and bind the result together 
			# into a matrix
		mts <- cbind(
				seasonal.stat.work(x[,1],
					fun=fun, win=win, giveStack=F,...),
				seasonal.stat.work(x[,2], 
					fun=fun, win=win, giveStack=F,...)
				)
		for(i in 2:(dim(x)[2]-1)){

			mts <- cbind(mts, 
				seasonal.stat.work(x[,i+1],
					fun=fun, win=win, giveStack=F,...))
				
		}
		mts <- as.data.frame(mts)
		names(mts) <- attributes(x)$dimnames[[2]]
		return(mts)
	}
	
}

seasonal.stat.work <- function(x,fun,win=frequency(x),giveStack=F,...){
	
	#calcualte a statistic for each period in a time series
	
	nyears <- length(unique(floor(time(x))))
	syear <- start(x)[1]

	times <- seq(syear,by=1/win,length.out=length(x))
	years <- floor(times)

	stacked <- matrix(NA,nrow=nyears,ncol=win)

	#Must use a loop because ts may not start at begining of year
	thisYear <- syear
	for(i in 1:nyears){

		thisRow <- x[years==thisYear]
		stacked[thisYear - syear + 1, 1:length(thisRow)] <- thisRow
		thisYear = thisYear + 1
	}

	if(giveStack)
		return(stacked)
	else
		return(as.vector(apply(stacked,2,fun,...)))
	
}

ts.annual.mean <- function(x){

	options(warn=-1)
	nyears <- length(unique(floor(time(x))))
	syear <- start(x)[1]
	nst = vector('numeric',nyears)

	for(i in 1:nyears)	nst[i] <- mean(window(x,syear - 1 + i,
										c(syear - 1 + i, frequency(x))))
   
	nst <- ts(nst, start = c(syear,1), frequency = 1) 
	return(nst)

}

wapply <- function(x, fun, win.len = length(x), ...){
	
	r <- (length(x) %% win.len)
	if(r > 0) x <- x[1:(length(x)-r)]
	stack <- matrix(x,nrow = win.len)
	return(apply(stack,2,fun,...))
	
}


ts.next.time <- function(x){
	en <- end(x)
	en <- if(en[2] == frequency(x))
			c(en[1] + 1, 1)
		else
			c(en[1], en[2] + 1)
	return(en)
}

ts.order.start <- function(x,start.per){

		f = frequency(x)
		ord <- if(start.per == 1) 
				start.per:f
		   else if(start.per == f) 
				c(f, 1:(f - 1))
		   else
				c(start.per:f, 1:(start.per - 1))
		return(ord)
}

myboxplot <- function (x, ..., range = 1.5, width = NULL, varwidth = FALSE, 
    notch = FALSE, outline = TRUE, names, boxwex = 0.8, plot = TRUE, 
    border = par("fg"), col = NULL, log = "", pars = NULL, horizontal = FALSE, 
    add = FALSE, at = NULL) 
{
    args <- list(x, ...)
    namedargs <- if (!is.null(attributes(args)$names)) 
        attributes(args)$names != ""
    else rep(FALSE, length.out = length(args))
    pars <- c(args[namedargs], pars)
    groups <- if (is.list(x)) 
        x
    else args[!namedargs]
    if (0 == (n <- length(groups))) 
        stop("invalid first argument")
    if (length(class(groups))) 
        groups <- unclass(groups)
    if (!missing(names)) 
        attr(groups, "names") <- names
    else {
        if (is.null(attr(groups, "names"))) 
            attr(groups, "names") <- 1:n
        names <- attr(groups, "names")
    }
    for (i in 1:n) groups[i] <- list(myboxplot.stats(groups[[i]], 
        range))
    stats <- matrix(0, nr = 5, nc = n)
    conf <- matrix(0, nr = 2, nc = n)
    ng <- out <- group <- numeric(0)
    ct <- 1
    for (i in groups) {
        stats[, ct] <- i$stats
        conf[, ct] <- i$conf
        ng <- c(ng, i$n)
        if ((lo <- length(i$out))) {
            out <- c(out, i$out)
            group <- c(group, rep.int(ct, lo))
        }
        ct <- ct + 1
    }
    z <- list(stats = stats, n = ng, conf = conf, out = out, 
        group = group, names = names)
    if (plot) {
        bxp(z, width, varwidth = varwidth, notch = notch, boxwex = boxwex, 
            border = border, col = col, log = log, pars = pars, 
            outline = outline, horizontal = horizontal, add = add, 
            at = at)
        invisible(z)
    }
    else z
}

#This function replaces boxplot.stats --- it uses the middle 50% and
#middle 90% instead of middle 50%(IQR) and 1.5*IQR.
myboxplot.stats <- function (x, coef = NULL, do.conf = TRUE, do.out =
TRUE)
{
  nna <- !is.na(x)
  n <- sum(nna)
  stats <- quantile(x, c(.05,.25,.5,.75,.95), na.rm = TRUE)
  iqr <- diff(stats[c(2, 4)])
  out <- x < stats[1] | x > stats[5]
  conf <- if (do.conf)
    stats[3] + c(-1.58, 1.58) * diff(stats[c(2, 4)])/sqrt(n)
  list(stats = stats, n = n, conf = conf, out = x[out & nna])
}

cfs2acft <- function(x){
	
	sts <- start(x)[1]
	sper <- start(x)[1]
	ets <- end(x)[1]
	eper <- end(x)[1]
	
	b <- 0 
	for(i in start(x)[1]:end(x)[1]){
		
		b <- b + 1
		r <- if(b == 1){
			sper:frequency(x)
		}else if(b == ets){
			 1:eper 
		}else{
			1:frequency(x)
		}
		
		for(j in r){
			
			
			
		}
			
	}
	
}


##########################
#
# End Functions
#
##########################