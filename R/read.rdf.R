read_rdf_header <- function(con,end){
	
	#assumes connection is already correctly positioned
	
	obj <- list()
	
	repeat{
		line <- readLines(con,n=1)
		if(line == end) break

		splitLine <- strsplit(line,':')[[1]]
		name <- splitLine[1]
		if(substr(splitLine[2],1,1) == ' ')
			splitLine[2] <- substr(splitLine[2],2,nchar(splitLine[2]))
		contents <- paste(splitLine[2:length(splitLine)],collapse=':')
		obj[[name]] <- contents
		
		if(end == 1) break
		
	}
	
	#returns the object, but connection is advanced
	return(obj)
	
}
	
read_rdf_meta <- function(rdf.con,rdf.obj){
	
	rdf.obj[['meta']] <- read_rdf_header(rdf.con,'END_PACKAGE_PREAMBLE')
	return(rdf.obj)
	
}

read_next_line_and_rewind <- function(con){
	
	str <- readLines(con,n=1)
	seek(con,seek(con)-nchar(str)-2)
	return(str)
	
}

read_rdf_run <- function(rdf.con,rdf.obj){
	
	this.run <- length(rdf.obj$runs) + 1
	rdf.obj$runs[[this.run]] <- read_rdf_header(rdf.con,'END_RUN_PREAMBLE')
	
		#time steps
	nts <- as.integer(rdf.obj$runs[[this.run]]$time_steps)
	rdf.obj$runs[[this.run]][['times']] <- readLines(rdf.con,n=nts)
	
		#Series
	nob <- 0
	repeat{
		
		nob <- nob + 1
		rdf.obj$runs[[this.run]][['objects']][[nob]] <- read_rdf_header(rdf.con,'END_SLOT_PREAMBLE')
		
			#name the objecst after their object.slot name
		obj.name <- rdf.obj$runs[[this.run]][['objects']][[nob]]$object_name
		slot.name <- rdf.obj$runs[[this.run]][['objects']][[nob]]$slot_name
		name <- paste(obj.name,slot.name,sep='.')
		names(rdf.obj$runs[[this.run]][['objects']])[nob] <- name
		
			#read in the extr two header pieces
		rdf.obj$runs[[this.run]][['objects']][[nob]]$units <- read_rdf_header(rdf.con,1)[[1]]
		rdf.obj$runs[[this.run]][['objects']][[nob]]$scale <- read_rdf_header(rdf.con,1)[[1]]

		rdf.obj$runs[[this.run]][['objects']][[nob]]$values <- as.numeric(readLines(rdf.con,n=nts))
		
			#END_COLUMN,END_SLOT, table slots need support here
		dummy <- readLines(rdf.con,n=2)
		
		if(read_next_line_and_rewind(rdf.con) == 'END_RUN'){
			dummy <- readLines(rdf.con,n=1)
			break
		}
	}
	return(rdf.obj)
	
}

read.rdf <- function(file){
	
	rdf.obj <- list()
	rdf.con <- file(file,'r')
	
	rdf.obj <- read_rdf_meta(rdf.con,rdf.obj)
	for(i in 1:as.numeric(rdf.obj$meta$number_of_runs))
		rdf.obj <- read_rdf_run(rdf.con,rdf.obj)
	
	close(rdf.con)
	
	return(rdf.obj)
}