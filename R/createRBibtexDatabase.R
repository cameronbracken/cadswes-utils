toBibFile <- function(pkgname,file,append=T){
	
	#writes a single bibtex package citation to a file
	
	if(append & !file.exists(file)) dummy <- file.create(file)
	
	entry <- toBibtex(citation(pkgname))
	entry <- if(is.list(entry)) as.character(entry[[1]]) else entry
	
	cat(entry,file=file,append=T,sep='\n')
	invisible()
}

createRBibtexDatabase <- function(file=NULL,append=T){
	# Writes a bibtex database file of all installed R packages
	# or returns a list of all citations
	#    if append=F an existing file will be overwritten
	
	installed <- as.list(rownames(installed.packages()))
	bib <- function(pkgname)toBibtex(citation(pkgname))
	
	if(is.null(file)){
		return(lapply(installed,bib))
	}else{
		if(!append ) dummy <- file.remove(file)
		dummy <- lapply(installed,toBibFile,file)
		invisible()
	}
}

#bibtexformat Library.bib -format -labels -f
