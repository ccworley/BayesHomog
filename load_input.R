# R Smiljanic Feb 2014 - Modified Sep 2016
# C Worley Aug 2018 - Modified for WG10
#def.par <- par(no.readonly = TRUE)
#
library(FITSio)
library(stringr)
library(gplots)
######################
#list.nodes <- c('Lumba','OACT','EPINARBO')  #,'IAC','MaxPlanck')
#this.release <- c('iDR6')
#node.template.file <- c('GES_iDR6_WG10_NodeTemplate_14062018_plus2ndExt.fits')
#node.files.path <- c('/Users/charlotteworley/Documents/GES/WG10/iDR6NodeQC/QCFiles/')
#
######################
#
# Function to load the results from the FITS files
#
load.nodes <- function(list.nodes,release=this.release,where.files=node.files.path) {
    if (length(list.nodes) == 0) stop('The list of nodes is empty')
  
    if (!is.character(list.nodes)) stop('The list of nodes should contain the Node names in character format')
#

    list.files <- paste0(where.files,'GES_',substr(release,1,4),'_WG10_',list.nodes,'.fits',sep="")
    # Use the first FITS file as reference for the numbers
    first.node <- readFITS(list.files[1], hdu=1, maxLines=25000)
    num.stars <- length(first.node$col[[1]])
    num.nodes <- length(list.nodes)
    #
    if (num.stars == 0) stop('There are no rows in the first FITS file')
    #
    # Copy the names of the files 
    #
    list.of.columns <- c('CNAME','GES_FLD','OBJECT','FILENAME','GES_TYPE','SETUP','SNR','VEL','E_VEL','VROT','E_VROT')
#    list.of.columns <- c('CNAME','GES_FLD','OBJECT','FILENAME','GES_TYPE','SETUP','SNR','VRAD','E_VRAD','VSINI','E_VSINI','REC_SETUP','REC_WG','WG')
    nodes.star.ids <- matrix(' ',num.stars,length(list.of.columns))
    colnames(nodes.star.ids) <- list.of.columns
    nodes.star.ids <- as.data.frame(nodes.star.ids)
    for (ik in seq(1,length(list.of.columns),1)) {
        numb.col <- which(first.node$colNames == list.of.columns[ik])
        nodes.star.ids[,ik] <- str_trim(first.node$col[[numb.col]]) # NODE.STAR.IDS is one and only comes from the FIRST NODE as reference
    }
    # Order by filename
    ii.files <- order(nodes.star.ids[,4])
    nodes.star.ids <- nodes.star.ids[ii.files,]
    # Now extract the parameters and errors
    list.param.names <- c('TEFF','E_TEFF','LOGG','E_LOGG','FEH','E_FEH','XI','E_XI','PECULI','REMARK','TECH')
    alt.list.names <- c('TEFF','E_TEFF','LOGG','E_LOGG','MH','E_MH','XI','E_XI','PECULI','REMARK','TECH')
    nodes.param <- array(NaN,c(num.stars,num.nodes,length(list.param.names)),dimnames=list(c(),list.nodes,list.param.names))
    for (ik in seq(1,num.nodes,1)) {
        node.to.extract <- readFITS(list.files[ik], hdu=1, maxLines=25000)
        col.files <- which(node.to.extract$colNames == "FILENAME")
        node.filenames <- str_trim(node.to.extract$col[[col.files]])
        
        if (num.stars != length(node.filenames)) stop(paste('The number of rows is different between ',list.files[1],' and ',list.files[ik],sep=""))

        # To guarantee that all Node results will be in the same order, and the order is by filename
        ii.order <- order(node.filenames)
        for (jk in seq(1,length(list.param.names),1)) {
            param.col <- which(node.to.extract$colNames == list.param.names[jk])
            if ((list.nodes[ik] == 'IAC')  && ((list.param.names[jk] == 'FEH') || (list.param.names[jk] == 'E_FEH'))) {
                param.col <- which(node.to.extract$colNames == alt.list.names[jk])
            }
            param.now <- node.to.extract$col[[param.col]]
            nodes.param[,ik,jk] <- param.now[ii.order]
        }
    }
    drops <- c('PECULI','REMARK','TECH')
    nodes.flags <- nodes.param[,,(dimnames(nodes.param)[[3]] %in% drops)]
    #
    # 'NaN' is not a flag
    #
    nodes.flags[nodes.flags == 'NaN'] <- " "
    #
    nodes.param <- nodes.param[,,(!dimnames(nodes.param)[[3]] %in% drops)]
    #
    for.output <- list(nodes.star.ids,nodes.param,nodes.flags)
    return(for.output)
}

#######################################
#
# Read the reference parameters of the benchmark stars
#
load.bench.fits <- function(fitsname,columns=columns.for.bench,path.file=path.for.bench,only.fgk=TRUE) {
    filename <- paste0(path.file,fitsname,sep="")
    bench.file <- readFITS(file=filename, hdu=1)
    list.of.columns <- columns
    num.stars <- length(bench.file$col[[1]])
    bench.stars <- matrix(' ',num.stars,length(list.of.columns))
    colnames(bench.stars) <- list.of.columns
    bench.stars <- as.data.frame(bench.stars)
    for (ik in seq(1,length(list.of.columns),1)) {
        numb.col <- which(bench.file$colNames == list.of.columns[ik])
        bench.stars[,ik] <- str_trim(bench.file$col[[numb.col]]) 
    }
    numeric.columns <- list.of.columns[!(list.of.columns %in% c('GES_FLD','GES_TYPE','GES_OBJECT','OBJECT','CNAME','FILENAME','CONSTFILES'))]
    for (each.col in numeric.columns) {
        ik <- which(list.of.columns == each.col)
        bench.stars[,ik] <- as.numeric(as.vector(bench.stars[,ik]))
    }
    if (only.fgk) {
        filter.bench.nans <- (is.nan(bench.stars$TEFF) | is.nan(bench.stars$LOGG) | (as.numeric(bench.stars$TEFF) > 9000) | (as.character(bench.stars$GES_FLD) == '32_Gem') | (as.character(bench.stars$GES_FLD) == 'HR1613') | (as.character(bench.stars$GES_FLD) == 'alf_Cep')) | (bench.stars$GES_TYPE != 'GE_SD_BM')
    } else {
        filter.bench.nans <- (is.nan(bench.stars$TEFF) | is.nan(bench.stars$LOGG))
    }
    bench.stars <- bench.stars[!filter.bench.nans,]
    return(bench.stars)
}

#par(def.par)#- reset to default

