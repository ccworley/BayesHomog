# R Smiljanic Feb 2014 - Modified Sep 2016 - Modified Jul, Aug, Sep, Oct 2018
# Plots for homog of iDR4 - Modified for iDR5 - and now for iDR6
#
######################
#
# Input
#
this.release <- c('iDR6')
list.columns.metadata <- c('CNAME','RA','DEC','GL','GB','JMAG','JMAGERR','HMAG','HMAGERR','KSMAG','KSMAGERR','UMAG', 'BMAG','VMAG','RMAG','IMAG',
                           'EBV','PMRAUCAC','PMRAERRUCAC','PMDECUCAC','PMDECERRUCAC')

list.columns.abun <- c("CNAME", "GESFIELD", "GESOBJECT", "CONSTFILES", "GESTYPE", "RECGRATINGS", "RECWG", "SNR", "VRAD", "VRADERR", "VSINI",
                       "VSINIERR", "TEFF", "TEFFERR", "LOGG", "LOGGERR", "FEH", "FEHERR", "XI", "XIERR", "NA1", "NA1ERR", "MG1", "MG1ERR", "AL1", "AL1ERR",
                       "SI1", "SI1ERR", "S1", "S1ERR", "CA1", "CA1ERR", "SC1", "SC1ERR", "SC2", "SC2ERR", "TI1", "TI1ERR", "TI2", "TI2ERR", "V1", "V1ERR",
                       "CR1", "CR1ERR", "CR2", "CR2ERR", "MN1", "MN1ERR", "FE1", "FE1ERR", "FE2", "FE2ERR", "CO1", "CO1ERR", "NI1", "NI1ERR", "CU1",
                       "CU1ERR", "ZN1", "ZN1ERR", "Y2", "Y2ERR", "ZR1", "ZR1ERR", "ZR2", "ZR2ERR", "BA2", "BA2ERR", "LA2", "LA2ERR",
                       "CE2", "CE2ERR", "ND2", "ND2ERR", "EU2", "EU2ERR")

#
######################
#
# Function to load the results from the FITS files
#
load.nodes <- function(list.nodes,release=this.release,where.files=node.files.path) {
  if (length(list.nodes) == 0) stop('The list of nodes is empty')
  if (!is.character(list.nodes)) stop('The list of nodes should contain the Node names in character format')
  #
  list.files <- paste0(where.files,'GES_',substr(release,1,4),'_WG10_',list.nodes,'.fits')
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
    if (list.nodes[ik] == 'Arcetri') {
      node.to.extract <- readFITS(list.files[ik], hdu=2, maxLines=25000)
    }
    col.files <- which(node.to.extract$colNames == "FILENAME")
    node.filenames <- str_trim(node.to.extract$col[[col.files]])
    if (list.nodes[ik] == 'ULB_OLD') {
      filter.uvl <- (substr(node.filenames,3,4) == "l_")
      node.filenames <- node.filenames[filter.uvl]
    }
    if (num.stars != length(node.filenames)) stop(paste('The number of rows is different between ',list.files[1],' and ',list.files[ik],sep=""))
    # To guarantee that all Node results will be in the same order, and the order is by filename
    ii.order <- order(node.filenames)
    for (jk in seq(1,length(list.param.names),1)) {
      param.col <- which(node.to.extract$colNames == list.param.names[jk])
      if (((list.nodes[ik] == 'IACb') || (list.nodes[ik] == 'IAC')) && ((list.param.names[jk] == 'FEH') || (list.param.names[jk] == 'E_FEH'))) {
        param.col <- which(node.to.extract$colNames == alt.list.names[jk])
      }
      if ((list.nodes[ik] == 'MyGIsFOS') && (list.param.names[jk] == 'FEH')) {
        param.col <- which(node.to.extract$colNames == 'FE1')
      }
      if ((list.nodes[ik] == 'MyGIsFOS') && (list.param.names[jk] == 'E_FEH')) {
        param.col <- which(node.to.extract$colNames == 'E_FE1')
      }
      param.now <- node.to.extract$col[[param.col]]
      if (list.nodes[ik] == 'ULB_OLD') {
        param.now <- param.now[filter.uvl]
      }
      if ((list.nodes[ik] == 'MyGIsFOS') && (list.param.names[jk] == 'FEH')) {
        param.now <- param.now-7.45
      }
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
# List of flags used by each node: Returns one list, where each component (for each node) is a list of 3 components (for each flag type)
#
list.flags <- function(nodes.flags) {
  list.of.nodes <- dimnames(nodes.flags)[[2]]
  list.of.flags <- dimnames(nodes.flags)[[3]]
  final.list.flags <- vector(mode="list",length=length(list.of.nodes))
  for (each.node in list.of.nodes) {
    node.list <- vector(mode="list",length=length(list.of.flags))
    col.node <- which(dimnames(nodes.flags)[[2]] == each.node)
    for (each.flag in list.of.flags) {
      col.flag <- which(dimnames(nodes.flags)[[3]] == each.flag)
      extracted.col <- nodes.flags[,col.node,col.flag]
      only.flags <- levels(factor(str_trim(extracted.col)))
      only.flags <- unlist(strsplit(only.flags,"|",fixed=TRUE))
      only.flags <- levels(factor(only.flags))
      only.flags <- unlist(strsplit(only.flags,"[ABC]"))
      only.flags <- levels(factor(substr(only.flags,1,11)))
      node.list[[col.flag]] <- only.flags
    }
    final.list.flags[[col.node]] <- node.list
  }
  names(final.list.flags) <- list.of.nodes
  for (i in seq(1,length(final.list.flags))) {
    names(final.list.flags[[i]]) <- list.of.flags
  }
  return(final.list.flags)
}
#######################################
#
# Functions to remove the flagged results and create the final recommended flags for iDR5
# flags.to.ignore.allnodes = stars with these flags, should be removed from ALL NODES
# flags.to.ignore.specific = stars with these flags, should have the Node specific result removed
#
# For iDR5 things get easier with the new dictionary, that has same sized flags always

eliminate.flagged.dr5 <- function(nodes.param,nodes.flags,flags.all,flags.specific) {
  #
  #print(str(nodes.param))
  list.nodes <- dimnames(nodes.flags)[[2]]
  filter.results <- vector("logical",length=dim(nodes.flags)[1])
  exclude.all <- vector("logical",length=dim(nodes.flags)[1])
  list.param <- dimnames(nodes.param)[[3]]
  col.peculi <- which(dimnames(nodes.flags)[[3]] == 'PECULI')
  col.tech <- which(dimnames(nodes.flags)[[3]] == 'TECH')
  
  #print('Number of NaN before flag eliminate')
  #print(sum(is.na(match(nodes.param[,,"TEFF"],"NaN"))))
  #print(sum(is.na(match(nodes.param[,'Lumba',"TEFF"],"NaN"))))
  #print(sum(is.na(match(nodes.param[,'EPINARBO',"TEFF"],"NaN"))))
  #print(sum(is.na(match(nodes.param[,'OACT',"TEFF"],"NaN"))))
  
  # Exclude in all nodes
  for (each.flag in flags.all) {
    #print(each.flag)
    for (each.node in list.nodes) {
      #print(each.node)
      node.col <- which(dimnames(nodes.flags)[[2]] == each.node)
      #print(node.col)
      # Only check PECULI and TECH
      exclude.all <- (exclude.all | grepl(each.flag,nodes.flags[,node.col,col.peculi]) | grepl(each.flag,nodes.flags[,node.col,col.tech]))
      #print(length(exclude.all))
      #print(sum(exclude.all))
      #print(which(exclude.all, arr.ind = FALSE, useNames = TRUE))
      #print('next')
      # exclude all will have TRUE in each place that was TRUE at least once
    }
  }
  # Exclude per node
  for (each.flag in flags.specific) {
    #print(each.flag)
    for (each.node in list.nodes) {
      node.col <- which(dimnames(nodes.flags)[[2]] == each.node)
      # Only check PECULI and TECH
      filter.results <- (exclude.all | grepl(each.flag,nodes.flags[,node.col,col.peculi]) | grepl(each.flag,nodes.flags[,node.col,col.tech]))
      #filter.results <- (grepl(each.flag,nodes.flags[,node.col,col.peculi]) | grepl(each.flag,nodes.flags[,node.col,col.tech]))
      # print(filter.results)
      # final filter.results will have TRUE in each place that was TRUE at least once + exclude.all
      node.col <- which(dimnames(nodes.param)[[2]] == each.node)
      for (each.param in list.param) {
        param.col <- which(dimnames(nodes.param)[[3]] == each.param)
        nodes.param[filter.results,node.col,param.col] <- "NaN"
      }
    }
  }
  #print('Number of NaN AFTER flag eliminate')
  #print(sum(is.na(match(nodes.param[,,"TEFF"],"NaN"))))
  #print(sum(is.na(match(nodes.param[,'Lumba',"TEFF"],"NaN"))))
  #print(sum(is.na(match(nodes.param[,'EPINARBO',"TEFF"],"NaN"))))
  #print(sum(is.na(match(nodes.param[,'OACT',"TEFF"],"NaN"))))
  
  #
  # Let's also produce the final combined flags for the recommended FITS file
  #
  num.stars <- length(nodes.flags[,1,1])
  num.flags <- length(nodes.flags[1,1,])
  recom.flags <- matrix("",num.stars,num.flags)
  list.flags <- dimnames(nodes.flags)[[3]]
  dimnames(recom.flags)[[2]] <- list.flags
  for (each.flag in list.flags) {
    flag.col <- which(dimnames(nodes.flags)[[3]] == each.flag)
    for (each.node in list.nodes) {
      node.col <- which(dimnames(nodes.flags)[[2]] == each.node)
      where.flags <- ((nchar(str_trim(nodes.flags[,node.col,flag.col])) > 1) & (!grepl("12099-11",str_trim(nodes.flags[,node.col,flag.col]))) & (!grepl("10399-11",str_trim(nodes.flags[,node.col,flag.col])))) # Flags "lack of time" are ignored
      recom.flags[where.flags,flag.col] <- str_trim(paste(recom.flags[where.flags,flag.col],nodes.flags[where.flags,node.col,flag.col],sep="|"))
    }
    where.first <- (substr(recom.flags[,flag.col],1,1) == "|")
    recom.flags[where.first,flag.col] <- substr(recom.flags[where.first,flag.col],2,nchar(recom.flags[where.first,flag.col]))
  }
  final.list <- list(nodes.param,recom.flags)
  return(final.list)
  #print(str(nodes.param))
  
}

recom.flags.per.cname.dr5 <- function(nodes.ids,recom.flags,final.abun.table) {
  # First, collapse all recom flags per cname
  # Then, check if there is a cname without param and/or abun that is without flags. If yes, list them and STOP    
  #
  list.flags <- dimnames(recom.flags)[[2]]
  for (each.cname in final.abun.table$CNAME) {
    which.recom.flags.rows <- which(nodes.ids$CNAME %in% each.cname)
    for (each.flag in list.flags) {
      flag.col.recom <- which(dimnames(recom.flags)[[2]] == each.flag)
      flag.col.abun <-  which(colnames(final.abun.table) == each.flag)
      new.composite.flag <- ""
      for (each.row in which.recom.flags.rows) {
        all.the.flags <- unlist(strsplit(str_trim(recom.flags[each.row,flag.col.recom]),"|",fixed=TRUE))
        #
        # Do not propagate some very specific flags
        #
        #list.of.specific.flags <- c('10309-11',
        #                            '10310-11','10311-11','10312-11','10313-11','10314-11','10315-11','10316-11','10317-11','10318-11',
        #                            '10319-11','10320-11','10390-11','10391-11','10500-11','13000-11','13002-11','13003-11','13010-11',
        #                            '13020-11','13021-11','13022-11','13023-11','13024-11','13025-11','13026-11','13027-11','13028-11',
        #                            '13029-11','13030-11','13031-11',paste0('1410',seq(1,9,1),'-11'),paste0('141',seq(10,99,1),'-11'),
        #                            paste0('1420',seq(1,9,1),'-11'),paste0('142',seq(10,99,1),'-11'),'10399-11','12099-11')
        #
        # New decision: just propagate everything!!!!!!!!!!!!!!!!!
        list.of.specific.flags <- c('0000000000000')
        #
        part.of.all.the.flags <- substr(all.the.flags,1,8)
        filter.the.specific <- (part.of.all.the.flags %in% list.of.specific.flags)
        all.the.flags <- all.the.flags[!filter.the.specific]
        if (length(all.the.flags) != 0) {
          for (each.indiv.flag in all.the.flags) {
            if (new.composite.flag == 'INDEF') { new.composite.flag <- "" }
            if ((nchar(each.indiv.flag) > 0) & !(each.indiv.flag %in% c('10301-11','10302-11','10303-11','10304-11','10305-11','10306-11','10307-11','10308-11'))) { # Leave some Node specific convergence flags without change
              this.flag <- each.indiv.flag #paste0(substr(each.indiv.flag,1,9),'00',substr(each.indiv.flag,12,16))# New decision, leave unchanged!!!!
              if (!grepl(this.flag,new.composite.flag)) {
                new.composite.flag <- paste(str_trim(new.composite.flag),str_trim(this.flag),sep="|")
              }
              if (substr(new.composite.flag,1,1) == "|") {
                new.composite.flag <- substr(new.composite.flag,2,nchar(new.composite.flag))
              }
            }
          }
        } else if (nchar(new.composite.flag) == 0) {
          new.composite.flag <- 'INDEF'
        }
      }
      # Find flags of SNR differing only by the SNR level, and adopt the least confident
      the.snr.flags <- c('10005','10010','10015','10020','10025','10030','10040','10050')
      if (nchar(new.composite.flag) > 16) {
        all.the.flags <- unlist(strsplit(str_trim(new.composite.flag),"|",fixed=TRUE))
        part.all.the.flags <- substr(all.the.flags,1,5)
        which.are.snr.flags <- (part.all.the.flags %in% the.snr.flags)
        sorted.snr.flags <- sort(all.the.flags[which.are.snr.flags])
        selected.snr.flag <- sorted.snr.flags[length(sorted.snr.flags)]
        all.the.flags[which.are.snr.flags] <- selected.snr.flag
        dup.flags <- duplicated(all.the.flags)
        new.composite.flag <- paste(all.the.flags[!dup.flags],collapse="|")
      }
      
      # Find repeated flags differing only by the confidence level, and adopt the least confident
      if (nchar(new.composite.flag) > 16) {
        all.the.flags <- unlist(strsplit(str_trim(new.composite.flag),"|",fixed=TRUE))
        all.flags.no.confidence.level <- substr(all.the.flags,1,14)
        all.flags.only.conf.level <- substr(all.the.flags,16,16)
        filter.dup.flags <- duplicated(all.flags.no.confidence.level)
        my.dup.flags <- all.flags.no.confidence.level[filter.dup.flags]
        for (each.case in my.dup.flags) {
          the.conf.levels <- all.flags.only.conf.level[all.flags.no.confidence.level == each.case]
          the.conf.levels <- sort(the.conf.levels)
          selected.conf.level <- the.conf.levels[length(the.conf.levels)]
          all.flags.only.conf.level[all.flags.no.confidence.level == each.case] <- selected.conf.level
        }
        #
        all.flags.without.repeat <- paste(all.flags.no.confidence.level[!filter.dup.flags],all.flags.only.conf.level[!filter.dup.flags],sep="-")
        all.flags.without.repeat <- sort(all.flags.without.repeat)
        #
        new.composite.flag <- paste0(all.flags.without.repeat,collapse="|")
      }
      if (max(nchar(recom.flags[which.recom.flags.rows,flag.col.recom])) > 0) {
        final.abun.table[(final.abun.table$CNAME == each.cname),flag.col.abun] <- new.composite.flag
      }
    }
  }
  # Blaim failure on lack of VROT
  filter.no.vrot <- (final.abun.table$VROT == '-999.90')
  final.abun.table$VROT[filter.no.vrot] <- 'INDEF'
  #    final.abun.table$TECH[filter.no.vrot] <- "10200-11-00-00-A"
  #
  # Blaim no abundances on lack of atmospheric parameters
  #
  filter.no.param <- (final.abun.table$TEFF == 'INDEF') & (final.abun.table$LOGG == 'INDEF') & (final.abun.table$FEH == 'INDEF') & (final.abun.table$XI == 'INDEF')
  #    final.abun.table$TECH[filter.no.param] <- "12003-11-00-05-A" # No abundance because no parameter was provided
  filter.no.logg <- (final.abun.table$TEFF != 'INDEF') & (final.abun.table$LOGG == 'INDEF')
  #    final.abun.table$TECH[filter.no.logg] <- "12003-11-00-02-A" # No abundance because no logg was provided
  #
  filter.no.feh <- (final.abun.table$TEFF != 'INDEF') & (final.abun.table$LOGG != 'INDEF') & (final.abun.table$FEH == 'INDEF')
  #
  #
  chem.elements <- colnames(final.abun.table[,seq(59,419,6)])
  num.elements <- length(chem.elements)
  filter.no.abun <- (apply((final.abun.table[,chem.elements] == "INDEF"),1,sum) == num.elements)
  #
  filter.no.abun.no.feh.no.flag <- filter.no.abun & filter.no.feh & (final.abun.table$TECH == 'INDEF') & (final.abun.table$PECULI == "INDEF") & (final.abun.table$REMARK == "INDEF")
  #    final.abun.table$TECH[filter.no.abun.no.feh.no.flag] <- "12003-11-00-03-A" # No abundance because no FeH was provided, but only add if there is no other flag there   
  #                                    
  filter.no.flag.no.abun <- (final.abun.table$PECULI == "INDEF") & (final.abun.table$TECH == 'INDEF') & (final.abun.table$REMARK == 'INDEF') & filter.no.abun
  #    if (sum(filter.no.flag.no.abun) > 0) {
  #        print('The following CNAMEs have no FLAGs and no Abundances:')
  #        print(final.abun.table$CNAME[filter.no.flag.no.abun])
  #        stop()
  #    }
  #
  
  return(final.abun.table)
}

######################################## 
#FIRST of all, correct IAC from the box grid results (teff > 7900 or < 3800 log g > 4.95 or < 0.05  feh > 0.8 or < -2.95)
#
correct.iacaip.grid <- function(in.nodes.param) {
  col.iac <- which(dimnames(in.nodes.param)[[2]] == "IACAIP")
  col.teff <- which(dimnames(in.nodes.param)[[3]] == "TEFF")
  col.logg <- which(dimnames(in.nodes.param)[[3]] == "LOGG")
  col.feh <- which(dimnames(in.nodes.param)[[3]] == "MH")
  selec.col <- in.nodes.param[,col.iac,col.teff]
  selec.col <- as.data.frame(selec.col)
  selec.col <- apply(selec.col,2,as.numeric)
  #
  filter.col1 <- !is.na(selec.col) & ((selec.col >= 7900) | (selec.col <= 3800))
  #
  selec.col <- in.nodes.param[,col.iac,col.logg]
  selec.col <- as.data.frame(selec.col)
  selec.col <- apply(selec.col,2,as.numeric)
  #    
  filter.col2 <- !is.na(selec.col) & ((selec.col <= 0.05) | (selec.col >= 4.95))
  #
  selec.col <- in.nodes.param[,col.iac,col.feh]
  selec.col <- as.data.frame(selec.col)
  selec.col <- apply(selec.col,2,as.numeric)
  #    
  filter.col3 <- !is.na(selec.col) & ((selec.col <= -2.95) | (selec.col >= 0.8))
  
  #
  final.filter <- filter.col1 | filter.col2 | filter.col3
  #
  for (each.param in dimnames(in.nodes.param)[[3]]) {
    col.param <- which(dimnames(in.nodes.param)[[3]] == each.param)
    in.nodes.param[final.filter,col.iac,each.param] <- "NaN"
  }
  
  
  
  return(in.nodes.param)
}

######################################## 
#FIRST of all, correct IAC, MaxPlanck from the box grid results (teff > 7900 or < 3800 log g > 4.95 or < 0.05  feh > 0.8 or < -2.95)
#
correct.iacmp.grid <- function(in.nodes.param) {
  col.iac <- which(dimnames(in.nodes.param)[[2]] == "IAC")
  col.teff <- which(dimnames(in.nodes.param)[[3]] == "TEFF")
  col.logg <- which(dimnames(in.nodes.param)[[3]] == "LOGG")
  col.feh <- which(dimnames(in.nodes.param)[[3]] == "MH")
  selec.col <- in.nodes.param[,col.iac,col.teff]
  selec.col <- as.data.frame(selec.col)
  selec.col <- apply(selec.col,2,as.numeric)
  #
  filter.col1 <- !is.na(selec.col) & ((selec.col >= 7900) | (selec.col <= 3800))
  #
  selec.col <- in.nodes.param[,col.iac,col.logg]
  selec.col <- as.data.frame(selec.col)
  selec.col <- apply(selec.col,2,as.numeric)
  #    
  filter.col2 <- !is.na(selec.col) & ((selec.col <= 0.05) | (selec.col >= 4.95))
  #
  selec.col <- in.nodes.param[,col.iac,col.feh]
  selec.col <- as.data.frame(selec.col)
  selec.col <- apply(selec.col,2,as.numeric)
  #    
  filter.col3 <- !is.na(selec.col) & ((selec.col <= -2.95) | (selec.col >= 0.8))
  
  #
  final.filter <- filter.col1 | filter.col2 | filter.col3
  #
  for (each.param in dimnames(in.nodes.param)[[3]]) {
    col.param <- which(dimnames(in.nodes.param)[[3]] == each.param)
    in.nodes.param[final.filter,col.iac,each.param] <- "NaN"
  }

  #Max Planck
  col.maxplanck <- which(dimnames(in.nodes.param)[[2]] == "MaxPlanck")
  col.teff <- which(dimnames(in.nodes.param)[[3]] == "TEFF")
  col.logg <- which(dimnames(in.nodes.param)[[3]] == "LOGG")
  col.feh <- which(dimnames(in.nodes.param)[[3]] == "FEH")
  selec.col <- in.nodes.param[,col.maxplanck,col.teff]
  selec.col <- as.data.frame(selec.col)
  selec.col <- apply(selec.col,2,as.numeric)
  #
  filter.col1 <- !is.na(selec.col) & ((selec.col >= 6950) | (selec.col <= 4000))
  #
  selec.col <- in.nodes.param[,col.maxplanck,col.logg]
  selec.col <- as.data.frame(selec.col)
  selec.col <- apply(selec.col,2,as.numeric)
  #    
  filter.col2 <- !is.na(selec.col) & ((selec.col <= 1.0) | (selec.col >= 4.90))
  #
  selec.col <- in.nodes.param[,col.maxplanck,col.feh]
  selec.col <- as.data.frame(selec.col)
  selec.col <- apply(selec.col,2,as.numeric)
  #    
  filter.col3 <- !is.na(selec.col) & ((selec.col <= -2.52) | (selec.col >= 0.6))
    #
  #
  selec.colt <- in.nodes.param[,col.maxplanck,col.teff]
  selec.colt <- as.data.frame(selec.colt)
  selec.colt <- apply(selec.colt,2,as.numeric)
  selec.coll <- in.nodes.param[,col.maxplanck,col.logg]
  selec.coll <- as.data.frame(selec.coll)
  selec.coll <- apply(selec.coll,2,as.numeric)
  
  filter.col4 <- !is.na(selec.colt) & !is.na(selec.coll) & selec.colt <= 4350 & selec.coll <= 4.0  & selec.coll >= 2.7
  #print(sum(filter.col4))
  #
  final.filter <- filter.col1 | filter.col2 | filter.col3 | filter.col4
  #print(sum(filter.col1 | filter.col2 | filter.col3))
  #print(sum(final.filter))
  #
  for (each.param in dimnames(in.nodes.param)[[3]]) {
    col.param <- which(dimnames(in.nodes.param)[[3]] == each.param)
    in.nodes.param[final.filter,col.maxplanck,each.param] <- "NaN"
  }
  
  #Max Planckb
  col.maxplanck <- which(dimnames(in.nodes.param)[[2]] == "MaxPlanckb")
  col.teff <- which(dimnames(in.nodes.param)[[3]] == "TEFF")
  col.logg <- which(dimnames(in.nodes.param)[[3]] == "LOGG")
  col.feh <- which(dimnames(in.nodes.param)[[3]] == "FEH")
  
  selec.col <- in.nodes.param[,col.maxplanck,col.teff]
  selec.col <- as.data.frame(selec.col)
  selec.col <- apply(selec.col,2,as.numeric)
  #
  filter.col1 <- !is.na(selec.col) & ((selec.col >= 6950) | (selec.col <= 4000))
  #
  selec.col <- in.nodes.param[,col.maxplanck,col.logg]
  selec.col <- as.data.frame(selec.col)
  selec.col <- apply(selec.col,2,as.numeric)
  #    
  filter.col2 <- !is.na(selec.col) & ((selec.col <= 1.0) | (selec.col >= 4.90))
  #
  selec.col <- in.nodes.param[,col.maxplanck,col.feh]
  selec.col <- as.data.frame(selec.col)
  selec.col <- apply(selec.col,2,as.numeric)
  #    
  filter.col3 <- !is.na(selec.col) & ((selec.col <= -2.52) | (selec.col >= 0.6))
  
  #
  selec.colt <- in.nodes.param[,col.maxplanck,col.teff]
  selec.colt <- as.data.frame(selec.colt)
  selec.colt <- apply(selec.colt,2,as.numeric)
  selec.coll <- in.nodes.param[,col.maxplanck,col.logg]
  selec.coll <- as.data.frame(selec.coll)
  selec.coll <- apply(selec.coll,2,as.numeric)
  
  filter.col4 <- !is.na(selec.colt) & !is.na(selec.coll) & selec.colt <= 4350 & selec.coll <= 4.0  & selec.coll >= 2.7
  #print(sum(filter.col4))
  #
  final.filter <- filter.col1 | filter.col2 | filter.col3 | filter.col4
  #print(sum(filter.col1 | filter.col2 | filter.col3))
  #print(sum(final.filter))
  #
  for (each.param in dimnames(in.nodes.param)[[3]]) {
    col.param <- which(dimnames(in.nodes.param)[[3]] == each.param)
    in.nodes.param[final.filter,col.maxplanck,each.param] <- "NaN"
  }
  
  return(in.nodes.param)
}


#######################################
#
# Find outlier results
#
# Right now flags maximum 1 Node per star, if at the same time Teff,logg AND Fe/H are different by more than 500 K or 0.5 dex in logg, 
# OR 0.5 dex in FEH from the MEAN

find.outliers <- function(in.nodes.param) {
  list.param <- c('TEFF','LOGG','FEH')
  limits.param <- c(500000,50,50)
  final.filter <- array(FALSE,c(dim(in.nodes.param)[1],dim(in.nodes.param)[2],3))
  for (each.param in list.param) {
    param.col <- which(dimnames(in.nodes.param[,,])[[3]] == each.param)
    param.as.numbers <- as.data.frame(in.nodes.param[,,param.col])
    param.as.numbers <- apply(param.as.numbers,2,as.numeric)
    mean.param <- apply(param.as.numbers,1,mean,na.rm=TRUE)
    sd.param <- apply(param.as.numbers,1,sd,na.rm=TRUE)
    diff.to.mean <- abs(param.as.numbers-mean.param)
    max.per.star <- apply(diff.to.mean,1,max,na.rm=TRUE)
    #
    filter.max.inf <- !is.infinite(max.per.star)
    filter.max <- (diff.to.mean == max.per.star)
    comb.filter.max <- filter.max & filter.max.inf
    #
    filter.diff <- (diff.to.mean >= limits.param[(param.col+1)/2])
    # sometimes there are many nodes that disagree with the mean (like for GJ205, 5 nodes disagreed with the mean, and the mean is wrong anyway...)        
    # so we will limit for the cases where there is only one clear Node-outlier
    sum.diff <- apply(filter.diff,1,sum,na.rm=TRUE)
    new.filter.diff <- (filter.diff) & (sum.diff == 1)
    #
    comb.filter <- (comb.filter.max & new.filter.diff)
    na.comb.filter <- is.na(comb.filter)
    final.filter[,,(param.col+1)/2] <- !na.comb.filter & comb.filter   
  }

  out.filter <- final.filter[,,1] | final.filter[,,2] | final.filter[,,3] 
  return(out.filter)
}

apply.filter.outliers <- function(in.nodes.param,filter.outliers) {
  list.nodes <- dimnames(in.nodes.param)[[2]]
  list.param <- dimnames(in.nodes.param)[[3]]
  for (each.node in list.nodes) {
    node.col <- which(dimnames(in.nodes.param[,,])[[2]] == each.node)
    filter.col <- filter.outliers[,node.col]
    ik <- 1
    while (ik <= dim(in.nodes.param)[[3]]) {
      in.nodes.param[filter.col,node.col,ik] <- "NaN"
      ik <- ik+1
    }
  }
  return(in.nodes.param)
}

#######################################
#
# Find outlier results compared to Accepted Values (WG11 xmat)
#
# Right now flags maximum 1 Node per star, if at the same time Teff,logg AND Fe/H are different by more than 500 K or 0.5 dex in logg, 
# OR 0.5 dex in FEH from the MEAN

find.outliers.accepted <- function(in.nodes.param,metadata.nodes.param,in.bench.param) {
  list.param <- c('TEFF','LOGG','FEH')
  limits.param <- c(50000,50,50)
  final.filter <- array(FALSE,c(dim(in.nodes.param)[1],dim(in.nodes.param)[2],3))
  for (each.param in list.param) {
    param.col <- which(dimnames(in.nodes.param[,,])[[3]] == each.param)
    param.as.numbers <- as.data.frame(in.nodes.param[,,param.col])
    param.as.numbers <- apply(param.as.numbers,2,as.numeric)
    
    if ("TEFF" %in% each.param) {
      bench.param.col <- in.bench.param$TEFF
    } else if ("LOGG" %in% each.param) {
      bench.param.col <- in.bench.param$LOGG
    } else if ("FEH" %in% each.param) {
      bench.param.col <- in.bench.param$FEH
    }
    
    # Find value of Accepted Star - match on CNAME
    mean.param <- vector('numeric',length=nrow(in.nodes.param))
    for (ik in seq(1,nrow(in.nodes.param))) {
      mean.param[ik] <- bench.param.col[which(in.bench.param$ID1 == metadata.nodes.param$CNAME[ik])]
      #print(which(in.bench.param$ID1 == metadata.nodes.param$CNAME[ik]))
      #print(bench.param.col[which(in.bench.param$ID1 == metadata.nodes.param$CNAME[ik])])
      #star.code[ik] <- which(bench.param$GES_FLD == metadata.of.bench.spectra$GES_FLD[ik])
    }
    
    #mean.param <- apply(param.as.numbers,1,mean,na.rm=TRUE)
    #sd.param <- apply(param.as.numbers,1,sd,na.rm=TRUE)
    diff.to.mean <- abs(param.as.numbers-mean.param)
    #print(diff.to.mean[350:700,])
    max.per.star <- apply(diff.to.mean,1,max,na.rm=TRUE)
    #
    filter.max.inf <- !is.infinite(max.per.star)
    #filter.max <- (diff.to.mean == max.per.star)
    comb.filter.max <- filter.max.inf #& filter.max 
    #
    filter.diff <- (diff.to.mean >= limits.param[(param.col+1)/2])
    # sometimes there are many nodes that disagree with the mean (like for GJ205, 5 nodes disagreed with the mean, and the mean is wrong anyway...)        
    # so we will limit for the cases where there is only one clear Node-outlier
    # sum.diff <- apply(filter.diff,1,sum,na.rm=TRUE)
    new.filter.diff <- (filter.diff) # & (sum.diff == 1)
    #
    comb.filter <- (comb.filter.max & new.filter.diff)
    na.comb.filter <- is.na(comb.filter)
    final.filter[,,(param.col+1)/2] <- !na.comb.filter & comb.filter   
  }
  
  out.filter <- final.filter[,,1] | final.filter[,,2] | final.filter[,,3] 
  return(out.filter)
}

#apply.filter.outliers <- functioFITn(in.nodes.param,filter.outliers) {
#  list.nodes <- dimnames(in.nodes.param)[[2]]
#  list.param <- dimnames(in.nodes.param)[[3]]
#  for (each.node in list.nodes) {
#    node.col <- which(dimnames(in.nodes.param[,,])[[2]] == each.node)
#    filter.col <- filter.outliers[,node.col]
#    ik <- 1
#    while (ik <= dim(in.nodes.param)[[3]]) {
#      in.nodes.param[filter.col,node.col,ik] <- "NaN"
#      ik <- ik+1
#    }
#  }
#  return(in.nodes.param)
#}


#######################################
#
# Read the reference parameters of the benchmark stars
#
load.benchwg11.fits <- function(fitsname,columns=columns.for.bench,path.file=path.for.bench,only.fgk=TRUE) {
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
  numeric.columns <- list.of.columns[!(list.of.columns %in% c('GES_FLD','GES_TYPE','ID1','GES_OBJECT','OBJECT','CNAME','FILENAME','CONSTFILES'))]
  for (each.col in numeric.columns) {
    ik <- which(list.of.columns == each.col)
    bench.stars[,ik] <- as.numeric(as.vector(bench.stars[,ik]))
  }
  if (only.fgk) {
    filter.bench.nans <- (is.nan(bench.stars$TEFF) | is.nan(bench.stars$LOGG) | (as.numeric(bench.stars$TEFF) > 8000) | (as.character(bench.stars$GES_FLD) == '32_Gem') | (as.character(bench.stars$GES_FLD) == 'HR1613') | (as.character(bench.stars$GES_FLD) == 'alf_Cep')) | (bench.stars$GES_TYPE != 'GE_SD_BM')
  } else {
    filter.bench.nans <- (is.nan(bench.stars$TEFF) | is.nan(bench.stars$LOGG) | (as.numeric(bench.stars$TEFF) > 8000) | (as.character(bench.stars$GES_FLD) == '32_Gem') | (as.character(bench.stars$GES_FLD) == 'HR1613') | (as.character(bench.stars$GES_FLD) == 'alf_Cep'))
  }
  bench.stars <- bench.stars[!filter.bench.nans,]
  return(bench.stars)
}
#######################################
#
# Plots of the Sun for the Nodes in iDR5
#
# parameters has to be in the format put out by load.nodes above: an array of (stars,nodes,atm.param); bench.param has to be in the format put out by load.bench above; 
# star.ids needs to be the full nodes.star.id from load.nodes.
#
sun.nodes <- function(parameters,star.ids,bench.param,release=substr(this.release,2,4)) {
  #    layout(cbind(c(1,2,3)))
  #
  #    Plot twice, separating 580 and 520 setups
  #    for (set.ik in levels(factor(wg11.nodes.id$SETUP))) {
  #    for (set.ik in c('U520','U580')) {
  for (set.ik in c('U520','U580')) {
    #
    filter.sun <- (((grepl('Sun',star.ids$OBJECT)) | (grepl('Sun',star.ids$OBJECT))) & (star.ids$SETUP == set.ik))
    if (length(dim(parameters)) == 3) {
      solar.param <- parameters[filter.sun,,]
      new.param <- array(NaN,dim(solar.param))
      dimnames(new.param) <- dimnames(solar.param)
    } else if (length(dim(parameters)) == 2) {
      if (sum(filter.sun) == 1) {
        solar.param <- parameters[filter.sun,]
        new.param <- vector("numeric",length(solar.param))
        names(new.param) <- names(solar.param)
      } else {
        solar.param <- parameters[filter.sun,]
        new.param <- matrix(NaN,dim(solar.param)[1],dim(solar.param)[2])
        dimnames(new.param) <- dimnames(solar.param)
      }
    } else {
      stop('Wrong number of dimensions in parameters file')
    }
    # New workaround to make it plot with clean.nodes.param, which is a character array
    
    if (sum(filter.sun) > 1) {
      if (!is.null(dim(solar.param))) {
        for (each.node in dimnames(solar.param)[[2]]) {
          node.col <- which(dimnames(solar.param)[[2]] == each.node)
          for (each.param in dimnames(solar.param)[[3]]) {
            param.col <- which(dimnames(solar.param)[[3]] == each.param)
            new.param[,node.col,param.col] <- apply((as.data.frame(solar.param[,node.col,param.col])),2,as.numeric)
          }
        }
      } else {
        stop('Have a look at names of solar.param')
      }
      solar.param <- new.param
    } else if (sum(filter.sun) == 1) {
      if (!is.null(dim(solar.param))) {
        for (each.param in dimnames(solar.param)[[2]]) {
          param.col <- which(dimnames(solar.param)[[2]] == each.param)
          new.param[,param.col] <- apply((as.data.frame(solar.param[,param.col])),2,as.numeric)
        }
      } else {
        for (each.param in names(solar.param)) {
          param.col <- which(names(solar.param) == each.param)
          new.param[param.col] <- apply((as.data.frame(solar.param[param.col])),2,as.numeric)
        }
      }
      solar.param <- new.param
    }
    #
    solar.ids <- star.ids[filter.sun,]
    filter.bench <- ((bench.param$GES_FLD== 'Sun') | (bench.param$GES_FLD== 'SUN'))
    ref.sun <- bench.param[filter.bench,]
    ref.sun$TEFF <- as.numeric(ref.sun$TEFF)
    ref.sun$E_TEFF <- as.numeric(ref.sun$E_TEFF)
    ref.sun$LOGG <- as.numeric(ref.sun$LOGG)
    ref.sun$E_LOGG<- as.numeric(ref.sun$E_LOGG)
    ref.sun$FEH <- as.numeric(ref.sun$FEH)
    ref.sun$XI <- as.numeric(ref.sun$XI)
    if (sum(filter.sun) > 1) {
      num.stars <- dim(solar.param)[1] # The number of different solar spectra
      num.nodes <- dim(solar.param)[2] # The number of Nodes
      node.names <- dimnames(solar.param)[[2]] # The Node names that will be copied to the x axis
    } else { # if there is only one star, then it is not an array anymore, but a simple matrix
      if (!is.null(dim(solar.param))) {
        num.stars <- sum(filter.sun) # The number of different solar spectra
        num.nodes <- dim(solar.param)[1] # The number of Nodes
        node.names <- dimnames(solar.param)[[1]] # The Node names that will be copied to the x axis
      } else {
        num.stars <- sum(filter.sun) # The number of different solar spectra
        num.nodes <- 1 # The number of Nodes
        node.names <- list.nodes[1] # The Node names that will be copied to the x axis
      }  
    }
    colors.to.use <- rainbow(sum(filter.sun)) # a vector of colors, one color for each solar spectrum
    #    
    x.plot <- seq(1,num.nodes,1) 
    #
    # Plot the Teff
    #
    out.plot <- paste('Plots/wg11_',release,'_',set.ik,'_sun_teff_nodes.eps',sep="")
    postscript(file=out.plot,horizontal=FALSE, onefile=FALSE, height=6, width=6, pointsize=10)
    par(oma=c(0.5,0.5,0.8,0.5))
    par(mar=c(5,6,2,2))
    par(mgp=c(4,1.2,0))
    #
    if (sum(filter.sun) > 1) {
      column.teff <- which(dimnames(solar.param)[[3]] == 'TEFF')
      teff.range <- range(solar.param[,,column.teff],na.rm=TRUE)
    } else {
      if (!is.null(dim(solar.param))) {
        column.teff <- which(dimnames(solar.param)[[2]] == 'TEFF')
        teff.range <- range(solar.param[,column.teff],na.rm=TRUE)
      } else {
        column.teff <- which(names(solar.param) == 'TEFF')
        teff.range <- range(solar.param[column.teff],na.rm=TRUE)
      }
    }
    if (min(teff.range,na.rm=TRUE) > ref.sun$TEFF-100) { teff.range[1] <- ref.sun$TEFF-100 }
    if (max(teff.range,na.rm=TRUE) < ref.sun$TEFF+100) { teff.range[2] <- ref.sun$TEFF+100 }
    #
    for (ik in seq(1,num.stars,1)) {
      if (sum(filter.sun) > 1) {
        y.plot <- solar.param[ik,,column.teff]
      } else {
        if (!is.null(dim(solar.param))) {
          y.plot <- solar.param[,column.teff]
        } else {
          y.plot <- solar.param[column.teff]
        }
      }
      if (ik == 1) {
        plotCI(x.plot,y.plot,pch=19,col=colors.to.use[ik],xlim=c(0,num.nodes+1),ylim=teff.range,xlab="",ylab=expression(paste('T'[eff],' (K)',sep="")), xaxt='n', yaxt='n',main=paste('WG11 - ',release,' - Sun - ',set.ik,sep=""), cex.main=1.2, cex.lab=1.5)
        legend("topleft",legend=paste('# ',sum(filter.sun),' spectra'),inset=0.02)
        box(which='plot',lwd=2)
        axis(side=1, at=seq(1,num.nodes,1), labels=node.names, lwd=2, las=2, col.axis='blue', cex.axis=0.9)
        #axis(side=1, at=seq(1,num.nodes,1), labels=FALSE, cex.axis=1.5, lwd=2)
        #
        if (diff(teff.range) < 500) {
          min.lim.teff <- (trunc(min(teff.range,na.rm=TRUE)/500)-1)*500
          max.lim.teff <- (ceiling(max(teff.range,na.rm=TRUE)/500)+1)*500
          step.lim.teff <- 100
        } else if (diff(teff.range) < 1000) {
          min.lim.teff <- (trunc(min(teff.range,na.rm=TRUE)/1000)-1)*1000
          max.lim.teff <- (ceiling(max(teff.range,na.rm=TRUE)/1000)+1)*1000
          step.lim.teff <- 200
        } else if (diff(teff.range) < 2500) {
          min.lim.teff <- (trunc(min(teff.range,na.rm=TRUE)/2500)-1)*2500
          max.lim.teff <- (ceiling(max(teff.range,na.rm=TRUE)/2500)+1)*2500
          step.lim.teff <- 500
        } else {
          stop(paste('The difference in Teff between multiple results is ',diff(teff.range)))
        }
        #
        axis(side=2, at=seq(min.lim.teff,max.lim.teff,step.lim.teff), labels=seq(min.lim.teff,max.lim.teff,step.lim.teff), cex.axis=1.5, lwd=2)
        axis(side=2, at=seq(min.lim.teff,max.lim.teff,step.lim.teff/10), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
        abline(h=ref.sun$TEFF,lty=1,col='red',lwd=2)
        text(0.3,ref.sun$TEFF-15,paste(round(ref.sun$TEFF,digits=0)),cex=1.2)
        abline(h=ref.sun$TEFF+100,lty=2,col='red',lwd=2)
        text(0.3,ref.sun$TEFF+85,paste('+100'),cex=1.2)
        abline(h=ref.sun$TEFF-100,lty=2,col='red',lwd=2)
        text(0.3,ref.sun$TEFF-115,paste('-100'),cex=1.2)
        #                     par(xpd=NA)
        #                text(x.plot,min(teff.range,na.rm=TRUE)-(4.5*step.lim.teff/10),node.names, col="blue", srt=90, cex=1.0)
        #                     par(xpd=FALSE)
      } else {
        if (sum(filter.sun) > 1) {
          plotCI(x.plot,y.plot,pch=19,col=colors.to.use[ik], xaxt='n', yaxt='n',add=TRUE)
        }
      }
    }
    if (sum(filter.sun) > 1) {
      for (jk in seq(1,num.nodes,1)) {
        boxplot(solar.param[,jk,column.teff], range=0, add=TRUE, at=jk, axes=FALSE, cex.lab=1.5, width=1.5)
      }
    }
    dev.off()
    #
    # Plot the logg
    #
    out.plot <- paste('Plots/wg11_',release,'_',set.ik,'_sun_logg_nodes.eps',sep="")
    postscript(file=out.plot,horizontal=FALSE, onefile=FALSE, height=6, width=6, pointsize=10)
    par(oma=c(0.5,0.5,0.8,0.5))
    par(mar=c(5,6,2,2))
    par(mgp=c(4,1.2,0))
    #
    if (sum(filter.sun) > 1) {
      column.logg <- which(dimnames(solar.param)[[3]] == 'LOGG')
      logg.range <- range(solar.param[,,column.logg],na.rm=TRUE)
    } else {
      if (!is.null(dim(solar.param))) {
        column.logg <- which(dimnames(solar.param)[[2]] == 'LOGG')
        logg.range <- range(solar.param[,column.logg],na.rm=TRUE)
      } else {
        column.logg <- which(names(solar.param) == 'LOGG')
        logg.range <- range(solar.param[column.logg],na.rm=TRUE)
      }
    }
    if (min(logg.range,na.rm=TRUE) > ref.sun$LOGG-0.2) { logg.range[1] <- ref.sun$LOGG-0.2 }
    if (max(logg.range,na.rm=TRUE) < ref.sun$LOGG+0.2) { logg.range[2] <- ref.sun$LOGG+0.2 }
    #
    for (ik in seq(1,num.stars,1)) {
      if (sum(filter.sun) > 1) {
        y.plot <- solar.param[ik,,column.logg]
      } else {
        if (!is.null(dim(solar.param))) {
          y.plot <- solar.param[,column.logg]
        } else {
          y.plot <- solar.param[column.logg]
        }
      }
      if (ik == 1) {
        plotCI(x.plot,y.plot,pch=19,col=colors.to.use[ik],xlim=c(0,num.nodes+1),ylim=logg.range,xlab="",ylab=expression(paste('log g (dex)',sep="")), xaxt='n', yaxt='n',main=paste('WG11 - ',release,' - Sun - ',set.ik,sep=""), cex.main=1.2, cex.lab=1.5)
        legend("topleft",legend=paste('# ',sum(filter.sun),' spectra'),inset=0.02)
        box(which='plot',lwd=2)
        axis(side=1, at=seq(1,num.nodes,1), labels=node.names, lwd=2, las=2, col.axis='blue', cex.axis=0.9)
        #            axis(side=1, at=seq(1,num.nodes,1), labels=FALSE, cex.axis=1.5, lwd=2)
        #
        if (diff(logg.range) < 0.5) {
          min.lim.logg <- (trunc(min(logg.range,na.rm=TRUE)/0.5)-1)*0.5
          max.lim.logg <- (ceiling(max(logg.range,na.rm=TRUE)/0.5)+1)*0.5
          step.lim.logg <- 0.1
        } else if (diff(logg.range) < 1.0) {
          min.lim.logg <- (trunc(min(logg.range,na.rm=TRUE)/1.0)-1)*1.0
          max.lim.logg <- (ceiling(max(logg.range,na.rm=TRUE)/1.0)+1)*1.0
          step.lim.logg <- 0.2
        } else if (diff(logg.range) < 2.5) {
          min.lim.logg <- (trunc(min(logg.range,na.rm=TRUE)/2.5)-1)*2.5
          max.lim.logg <- (ceiling(max(logg.range,na.rm=TRUE)/2.5)+1)*2.5
          step.lim.logg <- 0.5
        } else if (diff(logg.range) < 10) {
          min.lim.logg <- (trunc(min(logg.range,na.rm=TRUE)/10)-1)*10
          max.lim.logg <- (ceiling(max(logg.range,na.rm=TRUE)/10)+1)*10
          step.lim.logg <- 1.0
        } else {
          stop(paste('The difference in Logg between multiple results is ',diff(logg.range)))
        }
        #
        axis(side=2, at=seq(min.lim.logg,max.lim.logg,step.lim.logg), labels=seq(min.lim.logg,max.lim.logg,step.lim.logg), cex.axis=1.5, lwd=2)
        axis(side=2, at=seq(min.lim.logg,max.lim.logg,step.lim.logg/10), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
        abline(h=ref.sun$LOGG,lty=1,col='red',lwd=2)
        text(0.3,ref.sun$LOGG-0.02,paste(round(ref.sun$LOGG,digits=2)),cex=1.2)
        abline(h=ref.sun$LOGG+0.2,lty=2,col='red',lwd=2)
        text(0.3,ref.sun$LOGG+0.18,paste('+0.20'),cex=1.2)
        abline(h=ref.sun$LOGG-0.2,lty=2,col='red',lwd=2)
        text(0.3,ref.sun$LOGG-0.22,paste('-0.20'),cex=1.2)
        #                 par(xpd=NA)
        #            text(x.plot,min(logg.range,na.rm=TRUE)-(7*step.lim.logg/10),node.names, col="blue", srt=90, cex=1.0)
        #                 par(xpd=FALSE)
      } else {
        if (sum(filter.sun) > 1) {
          plotCI(x.plot,y.plot,pch=19,col=colors.to.use[ik], xaxt='n', yaxt='n',add=TRUE)
        }
      }
    }
    if (sum(filter.sun) > 1) {
      for (jk in seq(1,num.nodes,1)) {
        boxplot(solar.param[,jk,column.logg], range=0, add=TRUE, at=jk, axes=FALSE, cex.lab=1.5, width=1.5)
      }
    }
    dev.off()
    #
    # Plot the feh
    #
    out.plot <- paste('Plots/wg11_',release,'_',set.ik,'_sun_feh_nodes.eps',sep="")
    postscript(file=out.plot,horizontal=FALSE, onefile=FALSE, height=6, width=6, pointsize=10)
    par(oma=c(0.5,0.5,0.8,0.5))
    par(mar=c(5,6,2,2))
    par(mgp=c(4,1.2,0))
    #
    if (sum(filter.sun) > 1) {
      column.feh <- which(dimnames(solar.param)[[3]] == 'FEH')
      feh.range <- range(solar.param[,,column.feh],na.rm=TRUE)
    } else {
      if (!is.null(dim(solar.param))) {
        column.feh <- which(dimnames(solar.param)[[2]] == 'FEH')
        feh.range <- range(solar.param[,column.feh],na.rm=TRUE)
      } else {
        column.feh <- which(names(solar.param) == 'FEH')
        feh.range <- range(solar.param[column.feh],na.rm=TRUE)
      }
    }
    if (min(feh.range,na.rm=TRUE) > ref.sun$FEH-0.1) { feh.range[1] <- ref.sun$FEH-0.1 }
    if (max(feh.range,na.rm=TRUE) < ref.sun$FEH+0.1) { feh.range[2] <- ref.sun$FEH+0.1 }
    #
    for (ik in seq(1,num.stars,1)) {
      if (sum(filter.sun) > 1) {
        y.plot <- solar.param[ik,,column.feh]
      } else {
        if (!is.null(dim(solar.param))) {
          y.plot <- solar.param[,column.feh]
        } else {
          y.plot <- solar.param[column.feh]
        }
      }
      if (ik == 1) {
        plotCI(x.plot,y.plot,pch=19,col=colors.to.use[ik],xlim=c(0,num.nodes+1),ylim=feh.range,xlab="",ylab=expression(paste('[Fe/H] (dex)',sep="")), xaxt='n', yaxt='n',main=paste('WG11 - ',release,' - Sun - ',set.ik,sep=""), cex.main=1.2, cex.lab=1.5)
        legend("topleft",legend=paste('# ',sum(filter.sun),' spectra'),inset=0.02)
        box(which='plot',lwd=2)
        axis(side=1, at=seq(1,num.nodes,1), labels=node.names, lwd=2, las=2, col.axis='blue', cex.axis=0.9)
        #            axis(side=1, at=seq(1,num.nodes,1), labels=FALSE, cex.axis=1.5, lwd=2)
        #
        if (diff(feh.range) < 0.5) {
          min.lim.feh <- (trunc(min(feh.range,na.rm=TRUE)/0.5)-1)*0.5
          max.lim.feh <- (ceiling(max(feh.range,na.rm=TRUE)/0.5)+1)*0.5
          step.lim.feh <- 0.1
        } else if (diff(feh.range) < 1.0) {
          min.lim.feh <- (trunc(min(feh.range,na.rm=TRUE)/1.0)-1)*1.0
          max.lim.feh <- (ceiling(max(feh.range,na.rm=TRUE)/1.0)+1)*1.0
          step.lim.feh <- 0.2
        } else if (diff(feh.range) < 2.5) {
          min.lim.feh <- (trunc(min(feh.range,na.rm=TRUE)/2.5)-1)*2.5
          max.lim.feh <- (ceiling(max(feh.range,na.rm=TRUE)/2.5)+1)*2.5
          step.lim.feh <- 0.5
        } else {
          stop(paste('The difference in FEH between multiple results is ',diff(feh.range)))
        }
        #
        axis(side=2, at=seq(min.lim.feh,max.lim.feh,step.lim.feh), labels=seq(min.lim.feh,max.lim.feh,step.lim.feh), cex.axis=1.5, lwd=2)
        axis(side=2, at=seq(min.lim.feh,max.lim.feh,step.lim.feh/10), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
        abline(h=ref.sun$FEH,lty=1,col='red',lwd=2)
        text(0.3,ref.sun$FEH-0.01,paste(round(ref.sun$FEH,digits=2)),cex=1.2)
        abline(h=ref.sun$FEH+0.1,lty=2,col='red',lwd=2)
        text(0.3,ref.sun$FEH+0.08,paste('+0.10'),cex=1.2)
        abline(h=ref.sun$FEH-0.1,lty=2,col='red',lwd=2)
        text(0.3,ref.sun$FEH-0.12,paste('-0.10'),cex=1.2)
        #                 par(xpd=NA)
        #            text(x.plot,min(feh.range,na.rm=TRUE)-(6*step.lim.feh/10),node.names, col="blue", srt=90, cex=1.0)
        #                 par(xpd=FALSE)
      } else {
        if (sum(filter.sun) > 1) {
          plotCI(x.plot,y.plot,pch=19,col=colors.to.use[ik], xaxt='n', yaxt='n',add=TRUE)
        }
      }
    }
    if (sum(filter.sun) > 1) {
      for (jk in seq(1,num.nodes,1)) {
        boxplot(solar.param[,jk,column.feh], range=0, add=TRUE, at=jk, axes=FALSE, cex.lab=1.5, width=1.5)
      }
    }
    dev.off()
  }
}
#######################################
#
# Plots of the benchmark stars for the Nodes in iDRX
#
# parameters has to be in the format put out by load.nodes above: an array of (stars,nodes,atm.param); bench.param has to be in the format put out by load.bench above; 
# star.ids needs to be the full nodes.star.id from load.nodes.
#
# Here we look to match the CNAME and skip the Sun. The Sun has multiple cnames, so bet to do it in the function above
#
bench.nodes <- function(parameters,star.ids,bench.param,release=substr(this.release,2,4)) {
  #    layout(cbind(c(1,2,3)))
  #
  # New workaround to make it plot with clean.nodes.param, which is a character array
  new.param <- array(NaN,dim(parameters))
  dimnames(new.param) <- dimnames(parameters)
  if (length(dim(parameters)) == 3) {
    for (each.node in dimnames(parameters)[[2]]) {
      node.col <- which(dimnames(parameters)[[2]] == each.node)
      for (each.param in dimnames(parameters)[[3]]) {
        param.col <- which(dimnames(parameters)[[3]] == each.param)
        new.param[,node.col,param.col] <- apply((as.data.frame(parameters[,node.col,param.col])),2,as.numeric)
      }
    }
  } else {
    for (each.param in dimnames(parameters)[[2]]) {
      param.col <- which(dimnames(parameters)[[2]] == each.param)
      new.param[,param.col] <- as.numeric(parameters[,param.col])
    }
  }
  parameters <- new.param
  #
  bench.param$GES_FLD <- as.character(bench.param$GES_FLD)
  for (bench.id in bench.param$GES_FLD) {
    for (set.ik in c('U520','U580')) {
      if (!(bench.id %in% c('Sun','SUN','sun'))) {
        
        star.to.plot <- bench.param$GES_FLD[(bench.param$GES_FLD == bench.id)]
        filter.obj <- (grepl(star.to.plot,star.ids$GES_FLD) & (grepl(set.ik,star.ids$SETUP)))
        if (sum(filter.obj) != 0) {
          if (length(dim(parameters)) == 3) {
            from.bench.param <- parameters[filter.obj,,]
          } else {
            from.bench.param <- parameters[filter.obj,]
          }
          from.bench.ids <- star.ids[filter.obj,]
          filter.bench <- (bench.param$GES_FLD == bench.id)
          ref.from.bench <- bench.param[filter.bench,]
          ref.from.bench$TEFF <- as.numeric(ref.from.bench$TEFF)
          ref.from.bench$E_TEFF <- as.numeric(ref.from.bench$E_TEFF)
          ref.from.bench$LOGG <- as.numeric(ref.from.bench$LOGG)
          ref.from.bench$E_LOGG <- as.numeric(ref.from.bench$E_LOGG)
          ref.from.bench$FEH <- as.numeric(ref.from.bench$FEH)
          ref.from.bench$XI <- as.numeric(ref.from.bench$XI)
          #    
          if (sum(filter.obj) > 1) {
            num.stars <- dim(from.bench.param)[1] # The number of different from.bench spectra
            num.nodes <- dim(from.bench.param)[2] # The number of Nodes
            node.names <- dimnames(from.bench.param)[[2]] # The Node names that will be copied to the x axis
          } else if (sum(filter.obj) == 1) {
            num.stars <- 1 # The number of different from.bench spectra
            num.nodes <- dim(from.bench.param)[1] # The number of Nodes
            node.names <- dimnames(from.bench.param)[[1]] # The Node names that will be copied to the x axis
          } else {
            stop(paste('Something wrong, as no stars were found for ',bench.id))
          }
          #       
          colors.to.use <- rainbow(sum(filter.obj)) # a vector of colors, one color for each from.bench spectrum
          #
          x.plot <- seq(1,num.nodes,1) 
          #
          # Plot the Teff
          #
          if (sum(filter.obj) > 1) {
            column.teff <- which(dimnames(from.bench.param)[[3]] == 'TEFF')
            teff.range <- range(from.bench.param[,,column.teff],na.rm=TRUE)
          } else {
            column.teff <- which(dimnames(from.bench.param)[[2]] == 'TEFF')
            teff.range <- range(from.bench.param[,column.teff],na.rm=TRUE)
          }
          # Only plot in case there are values in teff.range, i.e., avoid if all values are missing
          if (!is.infinite(max(teff.range))) {
            
            output.plot <- paste('Plots/wg11_',release,'_',bench.id,'_',set.ik,'_teff_nodes.eps',sep="")
            postscript(output.plot,horizontal=FALSE, onefile=FALSE, height=6, width=6, pointsize=10)
            par(oma=c(0.5,0.5,0.8,0.5))
            par(mar=c(5,6,2,2))
            par(mgp=c(4,1.2,0))
            #
            #    
            #
            if (min(teff.range,na.rm=TRUE) > ref.from.bench$TEFF-100) { teff.range[1] <- ref.from.bench$TEFF-100 }
            if (max(teff.range,na.rm=TRUE) < ref.from.bench$TEFF+100) { teff.range[2] <- ref.from.bench$TEFF+100 }
            #
            for (ik in seq(1,num.stars,1)) {
              if (sum(filter.obj) > 1) {
                y.plot <- from.bench.param[ik,,column.teff]
              } else {
                y.plot <- from.bench.param[,column.teff]
              }
              if (ik == 1) {
                plotCI(x.plot,y.plot,pch=19,col=colors.to.use[ik],xlim=c(0,num.nodes+1),ylim=teff.range,xlab="",ylab=expression(paste('T'[eff],' (K)',sep="")), xaxt='n', yaxt='n',main=paste('WG11 - ',release,' - ',bench.id,' - ',set.ik,sep=""), cex.main=1.2, cex.lab=1.5)
                box(which='plot',lwd=2)
                legend("topleft",legend=paste('# ',sum(filter.obj),' spectra'),inset=0.02)
                axis(side=1, at=seq(1,num.nodes,1), labels=node.names, lwd=2, las=2, col.axis='blue', cex.axis=0.9)
                #            axis(side=1, at=seq(1,num.nodes,1), labels=FALSE, cex.axis=1.5, lwd=2)
                #
                if (diff(teff.range) < 500) {
                  min.lim.teff <- (trunc(min(teff.range,na.rm=TRUE)/500)-1)*500
                  max.lim.teff <- (ceiling(max(teff.range)/500)+1)*500
                  step.lim.teff <- 100
                } else if (diff(teff.range) < 1000) {
                  min.lim.teff <- (trunc(min(teff.range,na.rm=TRUE)/1000)-1)*1000
                  max.lim.teff <- (ceiling(max(teff.range,na.rm=TRUE)/1000)+1)*1000
                  step.lim.teff <- 200
                } else if (diff(teff.range) < 2500) {
                  min.lim.teff <- (trunc(min(teff.range,na.rm=TRUE)/2500)-1)*2500
                  max.lim.teff <- (ceiling(max(teff.range,na.rm=TRUE)/2500)+1)*2500
                  step.lim.teff <- 500
                } else if (diff(teff.range) < 5000) {
                  min.lim.teff <- (trunc(min(teff.range,na.rm=TRUE)/5000)-1)*5000
                  max.lim.teff <- (ceiling(max(teff.range,na.rm=TRUE)/5000)+1)*5000
                  step.lim.teff <- 1000
                } else {
                  stop(paste('The difference in Teff between multiple results is ',diff(teff.range)))
                }
                #
                axis(side=2, at=seq(min.lim.teff,max.lim.teff,step.lim.teff), labels=seq(min.lim.teff,max.lim.teff,step.lim.teff), cex.axis=1.5, lwd=2)
                axis(side=2, at=seq(min.lim.teff,max.lim.teff,step.lim.teff/10), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
                abline(h=ref.from.bench$TEFF,lty=1,col='red',lwd=2)
                text(0.3,ref.from.bench$TEFF-15,paste(round(ref.from.bench$TEFF,digits=0)),cex=1.2)
                abline(h=ref.from.bench$TEFF+100,lty=2,col='red',lwd=2)
                text(0.3,ref.from.bench$TEFF+85,paste('+100'),cex=1.2)
                abline(h=ref.from.bench$TEFF-100,lty=2,col='red',lwd=2)
                text(0.3,ref.from.bench$TEFF-115,paste('-100'),cex=1.2)
                #                 par(xpd=NA)
                #            text(x.plot,min(teff.range,na.rm=TRUE)-(5*step.lim.teff/10),node.names, col="blue", srt=90, cex=1.0)
                #                 par(xpd=FALSE)
              } else {
                plotCI(x.plot,y.plot,pch=19,col=colors.to.use[ik], xaxt='n', yaxt='n',add=TRUE)
              }
            }
            for (jk in seq(1,num.nodes,1)) {
              if (sum(filter.obj) > 1) {
                boxplot(from.bench.param[,jk,column.teff], range=0, add=TRUE, at=jk, axes=FALSE, cex.lab=1.5, width=1.5)
              } 
            }
            dev.off()
          }
          #
          # Plot the logg
          #
          if (sum(filter.obj) > 1) {
            column.logg <- which(dimnames(from.bench.param)[[3]] == 'LOGG')
            logg.range <- range(from.bench.param[,,column.logg],na.rm=TRUE)
          } else {
            column.logg <- which(dimnames(from.bench.param)[[2]] == 'LOGG')
            logg.range <- range(from.bench.param[,column.logg],na.rm=TRUE)
          }
          if (!is.infinite(max(logg.range))) {
            
            output.plot <- paste('Plots/wg11_',release,'_',bench.id,'_',set.ik,'_logg_nodes.eps',sep="")
            postscript(output.plot,horizontal=FALSE, onefile=FALSE, height=6, width=6, pointsize=10)
            par(oma=c(0.5,0.5,0.8,0.5))
            par(mar=c(5,6,2,2))
            par(mgp=c(4,1.2,0))
            #
            
            if (min(logg.range,na.rm=TRUE) > ref.from.bench$LOGG-0.2) { logg.range[1] <- ref.from.bench$LOGG-0.2 }
            if (max(logg.range,na.rm=TRUE) < ref.from.bench$LOGG+0.2) { logg.range[2] <- ref.from.bench$LOGG+0.2 }
            #
            for (ik in seq(1,num.stars,1)) {
              if (sum(filter.obj) > 1) {
                y.plot <- from.bench.param[ik,,column.logg]
              } else {
                y.plot <- from.bench.param[,column.logg]
              }
              if (ik == 1) {
                plotCI(x.plot,y.plot,pch=19,col=colors.to.use[ik],xlim=c(0,num.nodes+1),ylim=logg.range,xlab="",ylab=expression(paste('log g (dex)',sep="")), xaxt='n', yaxt='n',main=paste('WG11 - ',release,' - ',bench.id,' - ',set.ik,sep=""), cex.main=1.2, cex.lab=1.5)
                legend("topleft",legend=paste('# ',sum(filter.obj),' spectra'),inset=0.02)
                box(which='plot',lwd=2)
                axis(side=1, at=seq(1,num.nodes,1), labels=node.names, lwd=2, las=2, col.axis='blue', cex.axis=0.9)
                #            axis(side=1, at=seq(1,num.nodes,1), labels=FALSE, cex.axis=1.5, lwd=2)
                #
                if (diff(logg.range) < 0.5) {
                  min.lim.logg <- (trunc(min(logg.range,na.rm=TRUE)/0.5)-1)*0.5
                  max.lim.logg <- (ceiling(max(logg.range,na.rm=TRUE)/0.5)+1)*0.5
                  step.lim.logg <- 0.1
                } else if (diff(logg.range) < 1.0) {
                  min.lim.logg <- (trunc(min(logg.range,na.rm=TRUE)/1.0)-1)*1.0
                  max.lim.logg <- (ceiling(max(logg.range,na.rm=TRUE)/1.0)+1)*1.0
                  step.lim.logg <- 0.2
                } else if (diff(logg.range) < 2.5) {
                  min.lim.logg <- (trunc(min(logg.range,na.rm=TRUE)/2.5)-1)*2.5
                  max.lim.logg <- (ceiling(max(logg.range,na.rm=TRUE)/2.5)+1)*2.5
                  step.lim.logg <- 0.5
                } else if (diff(logg.range) < 5.0) {
                  min.lim.logg <- (trunc(min(logg.range,na.rm=TRUE)/5.0)-1)*5.0
                  max.lim.logg <- (ceiling(max(logg.range,na.rm=TRUE)/5.0)+1)*5.0
                  step.lim.logg <- 1.0
                } else if (diff(logg.range) < 10) {
                  min.lim.logg <- (trunc(min(logg.range,na.rm=TRUE)/10)-1)*10
                  max.lim.logg <- (ceiling(max(logg.range,na.rm=TRUE)/10)+1)*10
                  step.lim.logg <- 1.0
                } else {
                  stop(paste('The difference in Logg between multiple results is ',diff(logg.range)))
                }
                #
                axis(side=2, at=seq(min.lim.logg,max.lim.logg,step.lim.logg), labels=seq(min.lim.logg,max.lim.logg,step.lim.logg), cex.axis=1.5, lwd=2)
                axis(side=2, at=seq(min.lim.logg,max.lim.logg,step.lim.logg/10), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
                abline(h=ref.from.bench$LOGG,lty=1,col='red',lwd=2)
                text(0.3,ref.from.bench$LOGG-0.02,paste(round(ref.from.bench$LOGG,digits=2)),cex=1.2)
                abline(h=ref.from.bench$LOGG+0.2,lty=2,col='red',lwd=2)
                text(0.3,ref.from.bench$LOGG+0.18,paste('+0.20'),cex=1.2)
                abline(h=ref.from.bench$LOGG-0.2,lty=2,col='red',lwd=2)
                text(0.3,ref.from.bench$LOGG-0.22,paste('-0.20'),cex=1.2)
                #                 par(xpd=NA)
                #            text(x.plot,min(logg.range,na.rm=TRUE)-(7*step.lim.logg/10),node.names, col="blue", srt=90, cex=1.0)
                #                 par(xpd=FALSE)
              } else {
                plotCI(x.plot,y.plot,pch=19,col=colors.to.use[ik], xaxt='n', yaxt='n',add=TRUE)
              }
            }
            for (jk in seq(1,num.nodes,1)) {
              if (sum(filter.obj) > 1) {
                boxplot(from.bench.param[,jk,column.logg], range=0, add=TRUE, at=jk, axes=FALSE, cex.lab=1.5, width=1.5)
              }
            }
            dev.off()
          }
          #
          # Plot the feh
          #
          if (sum(filter.obj) > 1) {
            column.feh <- which(dimnames(from.bench.param)[[3]] == 'FEH')
            feh.range <- range(from.bench.param[,,column.feh],na.rm=TRUE)
          } else {
            column.feh <- which(dimnames(from.bench.param)[[2]] == 'FEH')
            feh.range <- range(from.bench.param[,column.feh],na.rm=TRUE)
          }
          if (!is.infinite(max(feh.range))) {
            
            output.plot <- paste('Plots/wg11_',release,'_',bench.id,'_',set.ik,'_feh_nodes.eps',sep="")
            postscript(output.plot,horizontal=FALSE, onefile=FALSE, height=6, width=6, pointsize=10)
            par(oma=c(0.5,0.5,0.8,0.5))
            par(mar=c(5,6,2,2))
            par(mgp=c(4,1.2,0))
            #
            if (!is.na(ref.from.bench$FEH)) {
              #        print(paste('This is feh.range = ',feh.range))
              if (min(feh.range,na.rm=TRUE) > ref.from.bench$FEH-0.1) { feh.range[1] <- ref.from.bench$FEH-0.1 }
              if (max(feh.range,na.rm=TRUE) < ref.from.bench$FEH+0.1) { feh.range[2] <- ref.from.bench$FEH+0.1 }
            }
            #
            for (ik in seq(1,num.stars,1)) {
              if (sum(filter.obj) > 1) {
                y.plot <- from.bench.param[ik,,column.feh]
              } else {
                y.plot <- from.bench.param[,column.feh]
              }        
              if (ik == 1) {
                plotCI(x.plot,y.plot,pch=19,col=colors.to.use[ik],xlim=c(0,num.nodes+1),ylim=feh.range,xlab="",ylab=expression(paste('[Fe/H] (dex)',sep="")), xaxt='n', yaxt='n',main=paste('WG11 - ',release,' - ',bench.id,' - ',set.ik,sep=""), cex.main=1.2, cex.lab=1.5)
                legend("topleft",legend=paste('# ',sum(filter.obj),' spectra'),inset=0.02)
                box(which='plot',lwd=2)
                axis(side=1, at=seq(1,num.nodes,1), labels=node.names, lwd=2, las=2, col.axis='blue', cex.axis=0.9)
                #            axis(side=1, at=seq(1,num.nodes,1), labels=FALSE, cex.axis=1.5, lwd=2)
                #
                if (diff(feh.range) < 0.5) {
                  min.lim.feh <- (trunc(min(feh.range,na.rm=TRUE)/0.5)-1)*0.5
                  max.lim.feh <- (ceiling(max(feh.range,na.rm=TRUE)/0.5)+1)*0.5
                  step.lim.feh <- 0.1
                } else if (diff(feh.range) < 1.0) {
                  min.lim.feh <- (trunc(min(feh.range,na.rm=TRUE)/1.0)-1)*1.0
                  max.lim.feh <- (ceiling(max(feh.range,na.rm=TRUE)/1.0)+1)*1.0
                  step.lim.feh <- 0.2
                } else if (diff(feh.range) < 2.5) {
                  min.lim.feh <- (trunc(min(feh.range,na.rm=TRUE)/2.5)-1)*2.5
                  max.lim.feh <- (ceiling(max(feh.range,na.rm=TRUE)/2.5)+1)*2.5
                  step.lim.feh <- 0.5
                } else if (diff(feh.range) < 5.0) {
                  min.lim.feh <- (trunc(min(feh.range,na.rm=TRUE)/5.0)-1)*5.0
                  max.lim.feh <- (ceiling(max(feh.range,na.rm=TRUE)/5.0)+1)*5.0
                  step.lim.feh <- 1.0
                } else {
                  stop(paste('The difference in FEH between multiple results is ',diff(feh.range)))
                }
                #
                axis(side=2, at=seq(min.lim.feh,max.lim.feh,step.lim.feh), labels=seq(min.lim.feh,max.lim.feh,step.lim.feh), cex.axis=1.5, lwd=2)
                axis(side=2, at=seq(min.lim.feh,max.lim.feh,step.lim.feh/10), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
                abline(h=ref.from.bench$FEH,lty=1,col='red',lwd=2)
                text(0.3,ref.from.bench$FEH-0.01,paste(round(ref.from.bench$FEH,digits=2)),cex=1.2)
                abline(h=ref.from.bench$FEH+0.1,lty=2,col='red',lwd=2)
                text(0.3,ref.from.bench$FEH+0.08,paste('+0.10'),cex=1.2)
                abline(h=ref.from.bench$FEH-0.1,lty=2,col='red',lwd=2)
                text(0.3,ref.from.bench$FEH-0.12,paste('-0.10'),cex=1.2)
                #                 par(xpd=NA)
                #            text(x.plot,min(feh.range,na.rm=TRUE)-(6*step.lim.feh/10),node.names, col="blue", srt=90, cex=1.0)
                #                 par(xpd=FALSE)
              } else {
                plotCI(x.plot,y.plot,pch=19,col=colors.to.use[ik], xaxt='n', yaxt='n',add=TRUE)
              }
            }
            for (jk in seq(1,num.nodes,1)) {
              if (sum(filter.obj) > 1) {
                boxplot(from.bench.param[,jk,column.feh], range=0, add=TRUE, at=jk, axes=FALSE, cex.lab=1.5, width=1.5)
              }
            }
            dev.off()
          }
        }
      }
    }
  }
}
######################################
#
# Plots of the benchmark stars for the Nodes in iDRX: Deltas Param 580 vs 520, Ref vs 580, Ref vs 520
#
compare.580and520.toref <- function(parameters,star.ids,bench.param,release=substr(this.release,2,4)) {
  # New workaround to make it plot with clean.nodes.param, which is a character array
  new.param <- array(NaN,dim(parameters))
  dimnames(new.param) <- dimnames(parameters)
  for (each.node in dimnames(parameters)[[2]]) {
    node.col <- which(dimnames(parameters)[[2]] == each.node)
    for (each.param in dimnames(parameters)[[3]]) {
      param.col <- which(dimnames(parameters)[[3]] == each.param)
      new.param[,node.col,param.col] <- apply((as.data.frame(parameters[,node.col,param.col])),2,as.numeric)
    }
  }
  parameters <- new.param
  #
  bench.param <- bench.param.to.numeric(bench.param)
  #
  bench.param$GES_FLD <- as.character(bench.param$GES_FLD)
  list.of.nodes <- dimnames(parameters)[[2]]
  for (each.node in list.of.nodes) {
    node.col <- which(dimnames(parameters)[[2]] == each.node)
    #        print(paste(each.node))
    plot.name <- paste('Plots/wg11_',release,'_bench_param_comparison_',each.node,'.eps',sep="")
    #        print(paste(plot.name))
    matrix.of.parameters <- matrix(NA,nrow=length(bench.param$GES_FLD),ncol=9) # ncol = 9 for Teff, logg, FEH, in Ref,U520 and U580
    colnames(matrix.of.parameters) <- c('Teff_Ref','logg_Ref','FeH_Ref','Teff_U520','logg_U520','FeH_U520','Teff_U580','logg_U580','FeH_U580')
    rownames(matrix.of.parameters) <- bench.param$GES_FLD
    matrix.of.parameters <- as.data.frame(matrix.of.parameters)
    for (each.star in bench.param$GES_FLD) {
      
      matrix.of.parameters$Teff_Ref[which(rownames(matrix.of.parameters) == each.star)] <- bench.param$TEFF[which(bench.param$GES_FLD == each.star)]
      matrix.of.parameters$logg_Ref[which(rownames(matrix.of.parameters) == each.star)] <- bench.param$LOGG[which(bench.param$GES_FLD == each.star)]
      matrix.of.parameters$FeH_Ref[which(rownames(matrix.of.parameters) == each.star)] <- bench.param$FEH[which(bench.param$GES_FLD == each.star)]
      filter.star.u520 <- (grepl(each.star,star.ids$GES_FLD)) & (star.ids$SETUP == 'U520')
      matrix.of.parameters$Teff_U520[which(rownames(matrix.of.parameters) == each.star)] <- mean(parameters[filter.star.u520,node.col,which(dimnames(parameters)[[3]] == 'TEFF')],na.rm=T)
      matrix.of.parameters$logg_U520[which(rownames(matrix.of.parameters) == each.star)] <- mean(parameters[filter.star.u520,node.col,which(dimnames(parameters)[[3]] == 'LOGG')],na.rm=T)
      matrix.of.parameters$FeH_U520[which(rownames(matrix.of.parameters) == each.star)] <- mean(parameters[filter.star.u520,node.col,which(dimnames(parameters)[[3]] == 'FEH')],na.rm=T)            
      filter.star.u580 <- (grepl(each.star,star.ids$GES_FLD)) & (star.ids$SETUP == 'U580')
      matrix.of.parameters$Teff_U580[which(rownames(matrix.of.parameters) == each.star)] <- mean(parameters[filter.star.u580,node.col,which(dimnames(parameters)[[3]] == 'TEFF')],na.rm=T)
      matrix.of.parameters$logg_U580[which(rownames(matrix.of.parameters) == each.star)] <- mean(parameters[filter.star.u580,node.col,which(dimnames(parameters)[[3]] == 'LOGG')],na.rm=T)
      matrix.of.parameters$FeH_U580[which(rownames(matrix.of.parameters) == each.star)] <- mean(parameters[filter.star.u580,node.col,which(dimnames(parameters)[[3]] == 'FEH')],na.rm=T)            
    }
    for (i in seq(1,ncol(matrix.of.parameters))) {
      matrix.of.parameters[(((matrix.of.parameters[,i]) == 'NaN') | is.na(matrix.of.parameters[,i])),i] <- NA
    }
    # Now do the plot
    postscript(file=plot.name,horizontal=FALSE, onefile=FALSE, height=6, width=12, pointsize=10)
    par(oma=c(0.5,0.5,0.8,0.5))
    par(mar=c(5,6,2,2))
    par(mgp=c(4,1.2,0))
    layout(rbind(c(1,2,3),c(4,5,6)))
    #
    # Teff 580
    #
    y.plot <- (matrix.of.parameters$Teff_U580)-(matrix.of.parameters$Teff_Ref)
    x.plot <- matrix.of.parameters$Teff_Ref
    x.range <- range(x.plot,na.rm=TRUE)
    x.range[1] <- x.range[1]-100
    x.range[2] <- x.range[2]+100
    y.range <- range(y.plot,na.rm=TRUE)
    y.range[1] <- y.range[1]-100
    y.range[2] <- y.range[2]+100
    if (is.infinite(x.range[1])) { x.range[[1]] <- 3000 }
    if (is.infinite(x.range[2])) { x.range[[2]] <- 7000 }
    if (is.infinite(y.range[1])) { y.range[[1]] <- -200 }
    if (is.infinite(y.range[2])) { y.range[[2]] <- 200 }
    #
    plotCI(x.plot,y.plot,pch=19,xlim=x.range,ylim=y.range,xlab=expression(paste('T'[eff],'_Ref (K)',sep="")),ylab=expression(paste('Delta(Ref-U580) (K)',sep="")), xaxt='n', yaxt='n', cex.main=1.2, cex.lab=1.5)
    box(which='plot',lwd=2)
    axis(side=1, at=seq(2000,9000,1000), labels=seq(2000,9000,1000), cex.axis=1.5, lwd=2)
    axis(side=1, at=seq(2000,9000,100), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    axis(side=2, at=seq(-9000,9000,200), labels=seq(-9000,9000,200), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(-9000,9000,50), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    abline(h=0,lwd=2,lty=2)
    x.line <- seq(2000,9000,10)
    if (sum(!is.na(y.plot)) > 3) {
      lin.reg <- lm(y.plot ~ x.plot)
      y.corr <- lin.reg$coefficients[1]+lin.reg$coefficients[2]*x.line
      lines(x.line,y.corr,lwd=2,lty=1,col='red')
      legend("topright", legend=bquote(paste(Delta, ' = (', .(round(lin.reg$coefficients[1], digits=2)) %+-% .(round(coef(summary(lin.reg))[1,"Std. Error"],digits=2)), ') + (',  .(round(lin.reg$coefficients[2], digits=2)) %+-% .(round(coef(summary(lin.reg))[2,"Std. Error"],digits=2)), ')*',T[eff],'(Ref)',sep="" )), lty=1, lwd=2, col='red', cex=0.7)
    }
    #
    # logg 580
    #
    y.plot <- (matrix.of.parameters$logg_U580)-(matrix.of.parameters$logg_Ref)
    x.plot <- matrix.of.parameters$logg_Ref
    x.range <- range(x.plot,na.rm=TRUE)
    x.range[1] <- x.range[1]-0.2
    x.range[2] <- x.range[2]+0.2
    y.range <- range(y.plot,na.rm=TRUE)
    y.range[1] <- y.range[1]-0.2
    y.range[2] <- y.range[2]+0.2
    if (is.infinite(x.range[1])) { x.range[[1]] <- 0.0 }
    if (is.infinite(x.range[2])) { x.range[[2]] <- 5.0 }
    if (is.infinite(y.range[1])) { y.range[[1]] <- -0.5 }
    if (is.infinite(y.range[2])) { y.range[[2]] <- 0.5 }
    #
    plotCI(x.plot,y.plot,pch=19,xlim=x.range,ylim=y.range,xlab=expression(paste('logg_Ref (dex)',sep="")),ylab=expression(paste('Delta(Ref-U580) (dex)',sep="")), xaxt='n', yaxt='n',main=paste('WG11 - ',release,' - ',each.node,sep=""), cex.main=1.2, cex.lab=1.5)
    box(which='plot',lwd=2)
    axis(side=1, at=seq(-2,7,1), labels=seq(-2,7,1), cex.axis=1.5, lwd=2)
    axis(side=1, at=seq(-2,7,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    axis(side=2, at=seq(-5,5,0.5), labels=seq(-5,5,0.5), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(-5,5,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    abline(h=0,lwd=2,lty=2)
    x.line <- seq(-2,7,0.1)
    if (sum(!is.na(y.plot)) > 3) {
      lin.reg <- lm(y.plot ~ x.plot)
      y.corr <- lin.reg$coefficients[1]+lin.reg$coefficients[2]*x.line
      lines(x.line,y.corr,lwd=2,lty=1,col='red')
      legend("topright", legend=bquote(paste(Delta, ' = (', .(round(lin.reg$coefficients[1], digits=2)) %+-% .(round(coef(summary(lin.reg))[1,"Std. Error"],digits=2)), ') + (',  .(round(lin.reg$coefficients[2], digits=2)) %+-% .(round(coef(summary(lin.reg))[2,"Std. Error"],digits=2)), ')*logg(Ref)',sep="" )), lty=1, lwd=2, col='red', cex=0.7)
    }
    #
    # FeH U580
    #
    y.plot <- (matrix.of.parameters$FeH_U580)-(matrix.of.parameters$FeH_Ref)
    x.plot <- matrix.of.parameters$FeH_Ref
    x.range <- range(x.plot,na.rm=TRUE)
    x.range[1] <- x.range[1]-0.1
    x.range[2] <- x.range[2]+0.1
    y.range <- range(y.plot,na.rm=TRUE)
    y.range[1] <- y.range[1]-0.1
    y.range[2] <- y.range[2]+0.1
    if (is.infinite(x.range[1])) { x.range[[1]] <- -3.0 }
    if (is.infinite(x.range[2])) { x.range[[2]] <- 0.5 }
    if (is.infinite(y.range[1])) { y.range[[1]] <- -0.5 }
    if (is.infinite(y.range[2])) { y.range[[2]] <- 0.5 }
    #
    plotCI(x.plot,y.plot,pch=19,xlim=x.range,ylim=y.range,xlab=expression(paste('[Fe/H]_Ref (dex)',sep="")),ylab=expression(paste('Delta(Ref-U580) (dex)',sep="")), xaxt='n', yaxt='n', cex.main=1.2, cex.lab=1.5)
    box(which='plot',lwd=2)
    axis(side=1, at=seq(-10,7,1), labels=seq(-10,7,1), cex.axis=1.5, lwd=2)
    axis(side=1, at=seq(-10,7,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    axis(side=2, at=seq(-5,4,0.2), labels=seq(-5,4,0.2), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(-5,4,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    abline(h=0,lwd=2,lty=2)
    x.line <- seq(-10,7,0.1)
    if (sum(!is.na(y.plot)) > 3) {
      lin.reg <- lm(y.plot ~ x.plot)
      y.corr <- lin.reg$coefficients[1]+lin.reg$coefficients[2]*x.line
      lines(x.line,y.corr,lwd=2,lty=1,col='red')
      legend("topright", legend=bquote(paste(Delta, ' = (', .(round(lin.reg$coefficients[1], digits=2)) %+-% .(round(coef(summary(lin.reg))[1,"Std. Error"],digits=2)), ') + (',  .(round(lin.reg$coefficients[2], digits=2)) %+-% .(round(coef(summary(lin.reg))[2,"Std. Error"],digits=2)), ')*[FeH](Ref)',sep="" )), lty=1, lwd=2, col='red', cex=0.7)
    }
    #
    #
    # Teff 520
    #
    y.plot <- (matrix.of.parameters$Teff_U520)-(matrix.of.parameters$Teff_Ref)
    x.plot <- matrix.of.parameters$Teff_Ref
    x.range <- range(x.plot,na.rm=TRUE)
    x.range[1] <- x.range[1]-100
    x.range[2] <- x.range[2]+100
    y.range <- range(y.plot,na.rm=TRUE)
    y.range[1] <- y.range[1]-100
    y.range[2] <- y.range[2]+100
    if (is.infinite(x.range[1])) { x.range[[1]] <- 3000 }
    if (is.infinite(x.range[2])) { x.range[[2]] <- 7000 }
    if (is.infinite(y.range[1])) { y.range[[1]] <- -200 }
    if (is.infinite(y.range[2])) { y.range[[2]] <- 200 }
    #
    plotCI(x.plot,y.plot,pch=19,xlim=x.range,ylim=y.range,xlab=expression(paste('T'[eff],'_Ref (K)',sep="")),ylab=expression(paste('Delta(Ref-U520) (K)',sep="")), xaxt='n', yaxt='n', cex.main=1.2, cex.lab=1.5)
    box(which='plot',lwd=2)
    axis(side=1, at=seq(2000,9000,1000), labels=seq(2000,9000,1000), cex.axis=1.5, lwd=2)
    axis(side=1, at=seq(2000,9000,100), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    axis(side=2, at=seq(-9000,9000,200), labels=seq(-9000,9000,200), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(-9000,9000,50), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    abline(h=0,lwd=2,lty=2)
    x.line <- seq(2000,9000,10)
    if (sum(!is.na(y.plot)) > 3) {
      lin.reg <- lm(y.plot ~ x.plot)
      y.corr <- lin.reg$coefficients[1]+lin.reg$coefficients[2]*x.line
      lines(x.line,y.corr,lwd=2,lty=1,col='red')
      legend("topright", legend=bquote(paste(Delta, ' = (', .(round(lin.reg$coefficients[1], digits=2)) %+-% .(round(coef(summary(lin.reg))[1,"Std. Error"],digits=2)), ') + (',  .(round(lin.reg$coefficients[2], digits=2)) %+-% .(round(coef(summary(lin.reg))[2,"Std. Error"],digits=2)), ')*',T[eff],'(Ref)',sep="" )), lty=1, lwd=2, col='red', cex=0.7)
    }
    #
    # logg 520
    #
    y.plot <- (matrix.of.parameters$logg_U520)-(matrix.of.parameters$logg_Ref)
    x.plot <- matrix.of.parameters$logg_Ref
    x.range <- range(x.plot,na.rm=TRUE)
    x.range[1] <- x.range[1]-0.2
    x.range[2] <- x.range[2]+0.2
    y.range <- range(y.plot,na.rm=TRUE)
    y.range[1] <- y.range[1]-0.2
    y.range[2] <- y.range[2]+0.2
    if (is.infinite(x.range[1])) { x.range[[1]] <- 0.0 }
    if (is.infinite(x.range[2])) { x.range[[2]] <- 5.0 }
    if (is.infinite(y.range[1])) { y.range[[1]] <- -0.5 }
    if (is.infinite(y.range[2])) { y.range[[2]] <- 0.5 }
    #        
    plotCI(x.plot,y.plot,pch=19,xlim=x.range,ylim=y.range,xlab=expression(paste('logg_Ref (dex)',sep="")),ylab=expression(paste('Delta(Ref-U520) (dex)',sep="")), xaxt='n', yaxt='n', cex.main=1.2, cex.lab=1.5)
    box(which='plot',lwd=2)
    axis(side=1, at=seq(-2,7,1), labels=seq(-2,7,1), cex.axis=1.5, lwd=2)
    axis(side=1, at=seq(-2,7,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    axis(side=2, at=seq(-5,5,0.5), labels=seq(-5,5,0.5), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(-5,5,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    abline(h=0,lwd=2,lty=2)
    x.line <- seq(-2,7,0.1)
    if (sum(!is.na(y.plot)) > 3) {
      lin.reg <- lm(y.plot ~ x.plot)
      y.corr <- lin.reg$coefficients[1]+lin.reg$coefficients[2]*x.line
      lines(x.line,y.corr,lwd=2,lty=1,col='red')
      legend("topright", legend=bquote(paste(Delta, ' = (', .(round(lin.reg$coefficients[1], digits=2)) %+-% .(round(coef(summary(lin.reg))[1,"Std. Error"],digits=2)), ') + (',  .(round(lin.reg$coefficients[2], digits=2)) %+-% .(round(coef(summary(lin.reg))[2,"Std. Error"],digits=2)), ')*logg(Ref)',sep="" )), lty=1, lwd=2, col='red', cex=0.7)
    }
    #
    # FeH U520
    #
    y.plot <- (matrix.of.parameters$FeH_U520)-(matrix.of.parameters$FeH_Ref)
    x.plot <- matrix.of.parameters$FeH_Ref
    x.range <- range(x.plot,na.rm=TRUE)
    x.range[1] <- x.range[1]-0.1
    x.range[2] <- x.range[2]+0.1
    y.range <- range(y.plot,na.rm=TRUE)
    y.range[1] <- y.range[1]-0.1
    y.range[2] <- y.range[2]+0.1
    if (is.infinite(x.range[1])) { x.range[[1]] <- -3.0 }
    if (is.infinite(x.range[2])) { x.range[[2]] <- 0.5 }
    if (is.infinite(y.range[1])) { y.range[[1]] <- -0.5 }
    if (is.infinite(y.range[2])) { y.range[[2]] <- 0.5 }
    #        
    plotCI(x.plot,y.plot,pch=19,xlim=x.range,ylim=y.range,xlab=expression(paste('[Fe/H]_Ref (dex)',sep="")),ylab=expression(paste('Delta(Ref-U520) (dex)',sep="")), xaxt='n', yaxt='n', cex.main=1.2, cex.lab=1.5)
    box(which='plot',lwd=2)
    axis(side=1, at=seq(-10,7,1), labels=seq(-10,7,1), cex.axis=1.5, lwd=2)
    axis(side=1, at=seq(-10,7,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    axis(side=2, at=seq(-5,4,0.2), labels=seq(-5,4,0.2), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(-5,4,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    abline(h=0,lwd=2,lty=2)
    x.line <- seq(-10,7,0.1)
    if (sum(!is.na(y.plot)) > 3) {
      lin.reg <- lm(y.plot ~ x.plot)
      y.corr <- lin.reg$coefficients[1]+lin.reg$coefficients[2]*x.line
      lines(x.line,y.corr,lwd=2,lty=1,col='red')
      legend("topright", legend=bquote(paste(Delta, ' = (', .(round(lin.reg$coefficients[1], digits=2)) %+-% .(round(coef(summary(lin.reg))[1,"Std. Error"],digits=2)), ') + (',  .(round(lin.reg$coefficients[2], digits=2)) %+-% .(round(coef(summary(lin.reg))[2,"Std. Error"],digits=2)), ')*[FeH](Ref)',sep="" )), lty=1, lwd=2, col='red', cex=0.7)
    }
    #
    dev.off()
  }
}
############################
#
# Plots of Teff vs logg per Node, divided by [Fe/H]
#
plot.teff.logg.nodes <- function(in.nodes.param,in.nodes.ids,release=substr(this.release,2,4)) {
  # Find out who are the GCs and OCs
  col.gestype <- which(dimnames(in.nodes.ids)[[2]] == "GES_TYPE")
  col.gesfld <- which(dimnames(in.nodes.ids)[[2]] == "GES_FLD")    
  types.ocs <- c("AR_CL","GE_CL","GE_SD_OC")
  types.gcs <- c("AR_SD_GC","GE_SD_GC")
  #
  filter.type <- (in.nodes.ids[,col.gestype] %in% types.ocs)
  names.ocs <- levels(factor(in.nodes.ids[filter.type,col.gesfld]))
  filter.type <- (in.nodes.ids[,col.gestype] %in% types.gcs)
  names.gcs <- levels(factor(in.nodes.ids[filter.type,col.gesfld]))
  # Correct M67
  names.gcs <- names.gcs[names.gcs != "M67"]
  if (sum(names.ocs == "M67") == 0) { names.ocs[length(names.ocs)+1] <- "M67" }
  # Did not use the GES TYPES directly, to be able to correct for M67
  filter.ocs <- in.nodes.ids[,col.gesfld] %in% names.ocs
  filter.gcs <- in.nodes.ids[,col.gesfld] %in% names.gcs
  filter.field <- !(in.nodes.ids[,col.gesfld] %in% names.ocs) & !(in.nodes.ids[,col.gesfld] %in% names.gcs)
  #
  plot.feh.limits <- function(max.feh,min.feh) {
    filter.feh <- (feh.vector <= max.feh) & (feh.vector > min.feh) & !is.na(feh.vector)
    use.gcs <- filter.gcs & filter.feh
    use.ocs <- filter.ocs & filter.feh
    use.field <- filter.field & filter.feh
    par(xpd=NA)
    plot(teff.vector[use.field],logg.vector[use.field], ylim=c(5,-0.5), xlim=c(7000,3300), ylab="log g (dex)", xlab=expression(paste(T[eff], " (K)")), pch="*", cex=2.5, lwd=3, cex.lab=1.7, cex.axis=1.7, axes=FALSE, xaxt='n', yaxt='n', col="black")
    text(5100,-0.2,paste(max.feh," >= [Fe/H] > ",min.feh), cex=1.7, col="black", lwd=1.7)
    par(xpd=FALSE)
    points(teff.vector[use.ocs],logg.vector[use.ocs], pch=3, col="red", cex=1.5, lwd=1)
    points(teff.vector[use.gcs],logg.vector[use.gcs], pch=1, col="blue", cex=1.5, lwd=2)
    str.max.feh <- 'none'
    if (max.feh == 0.5) { str.max.feh <- 'p0p50' }
    if (max.feh == -0.1) { str.max.feh <- 'm0p10' }
    if (max.feh == -0.7) { str.max.feh <- 'm0p70' }
    if (max.feh == -1.3) { str.max.feh <- 'm1p30' }
    if (max.feh == -1.9) { str.max.feh <- 'm1p90' }
    if (str.max.feh != 'none') {
      iso1 <- read.table(paste0('/work/rsmiljanic/Survey/DR5/Support/iso_parsec_12p5gyr_feh_',str.max.feh,'.dat'),header=F)
      iso2 <- read.table(paste0('/work/rsmiljanic/Survey/DR5/Support/iso_parsec_1p0gyr_feh_',str.max.feh,'.dat'),header=F)
      df.iso1 <- data.frame(x=10**(iso1$V6), y=iso1$V7)
      df.iso2 <- data.frame(x=10**(iso2$V6), y=iso2$V7)
      lines(df.iso1$x[iso1$V8 < 4],df.iso1$y[iso1$V8 < 4],col='orange',lwd=2,lty=1)
      lines(df.iso2$x[iso2$V8 < 4],df.iso2$y[iso2$V8 < 4],col='blueviolet',lwd=2,lty=1)
    }
    str.min.feh <- 'none'
    if (min.feh == -0.1) { str.min.feh <- 'm0p10' }
    if (min.feh == -0.7) { str.min.feh <- 'm0p70' }
    if (min.feh == -1.3) { str.min.feh <- 'm1p30' }
    if (min.feh == -1.9) { str.min.feh <- 'm1p90' }
    if (min.feh == -2.5) { str.min.feh <- 'm2p10' }
    if (str.min.feh != 'none') {
      iso1 <- read.table(paste0('/work/rsmiljanic/Survey/DR5/Support/iso_parsec_12p5gyr_feh_',str.min.feh,'.dat'),header=F)
      iso2 <- read.table(paste0('/work/rsmiljanic/Survey/DR5/Support/iso_parsec_1p0gyr_feh_',str.min.feh,'.dat'),header=F)
      df.iso1 <- data.frame(x=10**(iso1$V6), y=iso1$V7)
      df.iso2 <- data.frame(x=10**(iso2$V6), y=iso2$V7)
      lines(df.iso1$x[iso1$V8 < 4],df.iso1$y[iso1$V8 < 4],col='orange',lwd=2,lty=2)
      lines(df.iso2$x[iso2$V8 < 4],df.iso2$y[iso2$V8 < 4],col='blueviolet',lwd=2,lty=2)
    }
    axis(side=1, at=seq(7500,3000,-500), labels=seq(7500,3000,-500), cex.axis=1.5, lwd=2)
    axis(side=1, at=seq(7500,3000,-100), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    axis(side=2, at=seq(6,-1,-1), labels=seq(6,-1,-1), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(6,-1,-0.2), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    box(which="plot", lwd=2)
  }
  #
  if (length(dim(in.nodes.param)) == 3) {
    list.of.the.nodes <- dimnames(in.nodes.param)[[2]]
  } else {
    list.of.the.nodes <- c('Recommended_WG15')
  }
  for (each.node in list.of.the.nodes) {
    if (length(dim(in.nodes.param)) == 3) {
      col.node <- which(dimnames(in.nodes.param)[[2]] == each.node)
      out.plot <- paste('/work/rsmiljanic/Survey/DR6/Codes/Plots/wg11_',release,'_teff_logg_',each.node,'_by_feh.eps',sep="")
      teff.col <- which(dimnames(in.nodes.param)[[3]] == "TEFF")
      logg.col <- which(dimnames(in.nodes.param)[[3]] == "LOGG")
      feh.col <- which(dimnames(in.nodes.param)[[3]] == "FEH")
      teff.vector <- apply((as.data.frame(in.nodes.param[,col.node,teff.col])),2,as.numeric)
      logg.vector <- apply((as.data.frame(in.nodes.param[,col.node,logg.col])),2,as.numeric)
      feh.vector <- apply((as.data.frame(in.nodes.param[,col.node,feh.col])),2,as.numeric)
    } else {
      out.plot <- paste('/work/rsmiljanic/Survey/DR6/Codes/Plots/wg11_',release,'_teff_logg_',each.node,'_by_feh.eps',sep="")
      teff.col <- which(dimnames(in.nodes.param)[[2]] == "TEFF")
      logg.col <- which(dimnames(in.nodes.param)[[2]] == "LOGG")
      feh.col <- which(dimnames(in.nodes.param)[[2]] == "FEH")
      teff.vector <- as.numeric(in.nodes.param[,teff.col])
      logg.vector <- as.numeric(in.nodes.param[,logg.col])
      feh.vector <- as.numeric(in.nodes.param[,feh.col])
    }
    #
    postscript(file=out.plot,horizontal=FALSE, onefile=FALSE, height=6, width=12, pointsize=10)
    layout(rbind(c(1,2,3),c(4,5,6)))
    par(oma=c(1,2.5,1.0,0.5))
    par(mar=c(4,2.5,3,3))
    #  par(oma=c(0.5,0.5,0.8,0.5))
    #  par(mar=c(5,6,2,2))
    #  par(mgp=c(4,1.2,0))
    plot.title <- paste('UVES FGK-type: ',release,' ',each.node, sep=" ")
    #
    plot.feh.limits(0.50,-0.10)
    plot.feh.limits(-0.10,-0.70)
    par(xpd=NA)
    title(main=plot.title, cex.main=2.0)
    par(xpd=FALSE)
    plot.feh.limits(-0.70,-1.30)
    plot.feh.limits(-1.30,-1.90)
    plot.feh.limits(-1.90,-2.50)
    plot.feh.limits(-2.50,-3.10)
    #
    dev.off()
  }
}
############################
#
# Plots of Teff vs logg per Node, per cluster, with selection based on RV
#
plot.clusters.nodes <- function(in.nodes.param,in.nodes.ids,release=substr(this.release,2,4)) {
  # New workaround to make it plot with clean.nodes.param, which is a character array
  new.param <- array(NaN,dim(in.nodes.param))
  dimnames(new.param) <- dimnames(in.nodes.param)
  if (length(dim(in.nodes.param)) == 3) {
    for (each.node in dimnames(in.nodes.param)[[2]]) {
      node.col <- which(dimnames(in.nodes.param)[[2]] == each.node)
      for (each.param in dimnames(in.nodes.param)[[3]]) {
        param.col <- which(dimnames(in.nodes.param)[[3]] == each.param)
        new.param[,node.col,param.col] <- apply((as.data.frame(in.nodes.param[,node.col,param.col])),2,as.numeric)
      }
    }
  } else {
    for (each.param in dimnames(in.nodes.param)[[2]]) {
      param.col <- which(dimnames(in.nodes.param)[[2]] == each.param)
      new.param[,param.col] <- as.numeric(in.nodes.param[,param.col])
    }
  }
  in.nodes.param <- new.param
  #
  # Find out who are the GCs and OCs
  col.gestype <- which(dimnames(in.nodes.ids)[[2]] == "GES_TYPE")
  col.gesfld <- which(dimnames(in.nodes.ids)[[2]] == "GES_FLD")
  col.setup <- which(dimnames(in.nodes.ids)[[2]] == "SETUP")
  #
  col.vel <- which(dimnames(in.nodes.ids)[[2]] == "VEL")
  #  col.vel <- which(dimnames(in.nodes.ids)[[2]] == "VRAD")
  in.nodes.ids[,col.vel] <- apply((as.data.frame(in.nodes.ids[,col.vel])),2,as.numeric)
  #
  #types.clusters <- c("AR_CL","AR_SD_GC","AR_SD_OC","GE_CL","GE_SD_GC","GE_SD_OC")
  #
  #filter.type <- (in.nodes.ids[,col.gestype] %in% types.clusters)
  #names.clusters <- levels(factor(in.nodes.ids[filter.type,col.gesfld]))
  names.clusters <- c("Blanco1","Br25","Br31","Br36","Br44","Br81","gamma2_Vel","IC2391","IC2602","IC4665","M12","M15","M2","M67",
                      "Melotte71","NGC104","NGC1261","NGC1851","NGC1904","NGC2232","NGC2243","NGC2420","NGC2451","NGC2516","NGC2547",
                      "NGC2808","NGC3532","NGC362","NGC4372","NGC4815","NGC4833","NGC5927","NGC6005","NGC6253","NGC6259","NGC6553",
                      "NGC6633","NGC6705","NGC6752","NGC6802","Pismis18","Pleiades","Rup134","Trumpler20","Trumpler23")
  #
  if (length(dim(in.nodes.param)) == 3) {
    teff.col <- which(dimnames(in.nodes.param)[[3]] == "TEFF")
    logg.col <- which(dimnames(in.nodes.param)[[3]] == "LOGG")
  } else {
    teff.col <- which(dimnames(in.nodes.param)[[2]] == "TEFF")
    logg.col <- which(dimnames(in.nodes.param)[[2]] == "LOGG")
  }
  #
  plot.the.cluster <- function(cluster.name) {
    if (length(dim(in.nodes.param)) == 3) {
      x.plot <- in.nodes.param[filter.cluster,node.col,teff.col]
      y.plot <- in.nodes.param[filter.cluster,node.col,logg.col]
    } else {
      x.plot <- in.nodes.param[filter.cluster,teff.col]
      y.plot <- in.nodes.param[filter.cluster,logg.col]
    }
    rv.stars <- in.nodes.ids[filter.cluster,col.vel]
    if (sum(is.na(rv.stars)) != length(rv.stars)) {
      filter.members <- (rv.stars <= mean.rv+2*sd.rv) & (rv.stars >= mean.rv-2*sd.rv)
    } else {
      filter.members <- vector('logical',length=length(rv.stars))
      filter.members[!filter.members] <- T
    }
    x.members <- x.plot[filter.members]
    y.members <- y.plot[filter.members]
    x.non <- x.plot[!filter.members]
    y.non <- y.plot[!filter.members]
    if (sum(is.na(x.plot)) != length(x.plot)) {
      out.plot <- paste('/work/rsmiljanic/Survey/DR6/Codes/Plots/wg11_',release,'_iso_',each.cluster,'_',each.node,'_',each.setup,'.eps',sep="")
      iso.name <- paste("/work/rsmiljanic/Survey/DR5/Support/iso_",each.cluster,"_parsec.dat",sep="")
      postscript(file=out.plot,horizontal=FALSE, onefile=FALSE, height=6, width=6, pointsize=10)
      par(oma=c(0.5,0.5,0.8,0.5))
      par(mar=c(5,6,2,2))
      par(mgp=c(4,1.2,0))
      range.teff <- range(x.plot,na.rm=TRUE)
      range.teff <- c(range.teff[2],range.teff[1])
      range.teff[1] <- range.teff[1]+1000
      range.teff[2] <- range.teff[2]-500
      range.logg <- range(y.plot,na.rm=TRUE)
      range.logg <- c(range.logg[2],range.logg[1])
      range.logg[1] <- range.logg[1]+1.0
      range.logg[2] <- range.logg[2]-1.0
      plotCI(x.members,y.members, ylim=range.logg, xlim=range.teff, ylab="log g (dex)", xlab=expression(paste(T[eff], " (K)")), pch=1, cex=2.5, lwd=3, cex.lab=1.7, cex.axis=1.7, axes=FALSE, xaxt='n', yaxt='n', col="blue", main=paste('WG11 - ',release,' - ',each.cluster,' - ',each.node,sep=""))
      box(lwd=2)
      if (length(x.non) != 0) {
        points(x.non,y.non, pch=3, col="red", cex=1.5, lwd=2)
        num.stars <- sum(!is.na(x.members))
        num.non <- sum(!is.na(x.non))
        leg.text <- c(paste(num.stars,' (of ',length(x.members),' member stars)',sep=""),paste(num.non,' (of ',length(x.non),'non-member stars)',sep=""))
        legend("topleft",legend=leg.text,inset=0.02,pch=c(1,3),col=c('blue','red'))
      } else {
        num.stars <- sum(!is.na(x.members))
        leg.text <- c(paste(num.stars,' (of ',length(x.members),' member stars)',sep=""))
        legend("topleft",legend=leg.text,inset=0.02,pch=1,col='blue')
      }
      axis(side=1, at=seq(10000,2000,-500), labels=seq(10000,2000,-500), cex.axis=1.5, lwd=2)
      axis(side=1, at=seq(10000,2000,-100), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
      axis(side=2, at=seq(10.0,-2.0,-0.5), labels=seq(10.0,-2.0,-0.5), cex.axis=1.5, lwd=2)
      axis(side=2, at=seq(10.0,-2.0,-0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
      iso.one <- read.table(iso.name)
      teff.one <- 10**(iso.one$V6)
      lines(teff.one, iso.one$V7, lwd=1.2, lty=2)
      dev.off()
    }
  }
  #
  for (each.cluster in names.clusters) {
    for (each.setup in c("U520","U580")) {
      filter.cluster <- ((in.nodes.ids[,col.gesfld] == each.cluster) & (grepl(each.setup,in.nodes.ids[,col.setup])))
      if (sum(filter.cluster) != 0) { # There are stars in that cluster observed with that setup
        if (sum(is.na(in.nodes.ids[filter.cluster,col.vel])) != length(in.nodes.ids[filter.cluster,col.vel])) { # There are some RVs there
          mean.rv <- mean(in.nodes.ids[filter.cluster,col.vel],na.rm=TRUE)
          sd.rv <- sd(in.nodes.ids[filter.cluster,col.vel],na.rm=TRUE)
        } else { # Else, everybody is a member
          mean.rv = 0.0
          sd.rv = 10000.0
        }
        if (length(dim(in.nodes.param)) == 3) {
          for (each.node in dimnames(in.nodes.param)[[2]]) {
            node.col <- which(dimnames(in.nodes.param)[[2]] == each.node)
            plot.the.cluster(each.cluster)
          }
        } else {
          for (each.node in c('Recommended_WG11')) {
            plot.the.cluster(each.cluster)
          }
        }
      }
    }
  }
}
#
############################
#
# Plots of correlations between Fe/H and Teff or logg per Node, per cluster, with selection based on RV
#
plot.metal.clusters.nodes <- function(in.nodes.param,in.nodes.ids,release=substr(this.release,2,4)) {
  # New workaround to make it plot with clean.nodes.param, which is a character array
  new.param <- array(NaN,dim(in.nodes.param))
  dimnames(new.param) <- dimnames(in.nodes.param)
  if (length(dim(in.nodes.param)) == 3) {
    for (each.node in dimnames(in.nodes.param)[[2]]) {
      node.col <- which(dimnames(in.nodes.param)[[2]] == each.node)
      for (each.param in dimnames(in.nodes.param)[[3]]) {
        param.col <- which(dimnames(in.nodes.param)[[3]] == each.param)
        new.param[,node.col,param.col] <- apply((as.data.frame(in.nodes.param[,node.col,param.col])),2,as.numeric)
      }
    }
  } else {
    for (each.param in dimnames(in.nodes.param)[[2]]) {
      param.col <- which(dimnames(in.nodes.param)[[2]] == each.param)
      new.param[,param.col] <- as.numeric(in.nodes.param[,param.col])
    }
  }
  in.nodes.param <- new.param
  #
  
  
  
  #
  # Find out who are the GCs and OCs
  col.gestype <- which(dimnames(in.nodes.ids)[[2]] == "GES_TYPE")
  col.gesfld <- which(dimnames(in.nodes.ids)[[2]] == "GES_FLD")
  col.setup <- which(dimnames(in.nodes.ids)[[2]] == "SETUP")
  #
  col.vel <- which(dimnames(in.nodes.ids)[[2]] == "VEL")
  in.nodes.ids[,col.vel] <- apply((as.data.frame(in.nodes.ids[,col.vel])),2,as.numeric)
  #
  types.clusters <- c("AR_CL","GE_CL","GE_SD_OC","AR_SD_GC","GE_SD_GC")
  #
  filter.type <- (in.nodes.ids[,col.gestype] %in% types.clusters)
  names.clusters <- levels(factor(in.nodes.ids[filter.type,col.gesfld]))
  #
  #
  if (length(dim(in.nodes.param)) == 3) {
    teff.col <- which(dimnames(in.nodes.param)[[3]] == "TEFF")
    logg.col <- which(dimnames(in.nodes.param)[[3]] == "LOGG")
    feh.col <- which(dimnames(in.nodes.param)[[3]] == "FEH")
  } else {
    teff.col <- which(dimnames(in.nodes.param)[[2]] == "TEFF")
    logg.col <- which(dimnames(in.nodes.param)[[2]] == "LOGG")
    feh.col <- which(dimnames(in.nodes.param)[[2]] == "FEH")
  }
  #
  #
  plot.the.cluster <- function(cluster.name) {
    if (length(dim(in.nodes.param)) == 3) {
      teff.plot <- in.nodes.param[filter.cluster,node.col,teff.col]
      logg.plot <- in.nodes.param[filter.cluster,node.col,logg.col]
      feh.plot <- in.nodes.param[filter.cluster,node.col,feh.col]
    } else {
      teff.plot <- in.nodes.param[filter.cluster,teff.col]
      logg.plot <- in.nodes.param[filter.cluster,logg.col]
      feh.plot <- in.nodes.param[filter.cluster,feh.col]
    }
    rv.stars <- in.nodes.ids[filter.cluster,col.vel]
    filter.members <- (rv.stars <= mean.rv+2*sd.rv) & (rv.stars >= mean.rv-2*sd.rv)
    teff.members <- teff.plot[filter.members]
    logg.members <- logg.plot[filter.members]
    feh.members <- feh.plot[filter.members]
    teff.non <- teff.plot[!filter.members]
    logg.non <- logg.plot[!filter.members]
    feh.non <- feh.plot[!filter.members]
    if (sum(is.na(teff.plot)) != length(teff.plot)) {
      out.plot <- paste('/work/rsmiljanic/Survey/DR6/Codes/Plots/wg11_',release,'_metal_',each.cluster,'_',each.node,'_',each.setup,'.eps',sep="")
      postscript(file=out.plot,horizontal=FALSE, onefile=FALSE, height=6, width=18, pointsize=10)
      layout(rbind(c(1,2)))
      par(oma=c(1.0,0.5,1.0,0.5))
      par(mar=c(5,6,2,2))
      #            par(mgp=c(4,1.2,0))
      range.teff <- range(teff.plot,na.rm=TRUE)
      range.teff[1] <- range.teff[1]-200
      range.teff[2] <- range.teff[2]+200
      range.logg <- range(logg.plot,na.rm=TRUE)
      range.logg[1] <- range.logg[1]-0.2
      range.logg[2] <- range.logg[2]+0.2
      range.feh <- range(feh.plot,na.rm=TRUE)
      range.feh[1] <- range.feh[1]-0.2
      range.feh[2] <- range.feh[2]+0.2
      if ((range.feh[2]-range.feh[1]) < 1.0) {
        range.feh[1] <- range.feh[1]-0.3
        range.feh[2] <- range.feh[2]+0.3
      }
      #
      plotCI(teff.members,feh.members, ylim=range.feh, xlim=range.teff, ylab="[Fe/H] (dex)", xlab=expression(paste(T[eff], " (K)")), pch=1, cex=2.5, lwd=3, cex.lab=1.7, cex.axis=1.7, cex.main=1.5, axes=FALSE, xaxt='n', yaxt='n', col="blue")
      box(lwd=2)
      par(xpd=NA)
      title(main=paste('WG11 - ',release,' - ',each.cluster,' - ',each.node,sep=""),outer=TRUE)
      par(xpd=FALSE)
      if (length(teff.non) != 0) {
        points(teff.non,feh.non, pch=3, col="red", cex=1.5, lwd=2)
        num.stars <- sum(!is.na(teff.members))
        num.non <- sum(!is.na(teff.non))
        leg.text <- c(paste(num.stars,' (of ',length(teff.members),' member stars)',sep=""),paste(num.non,' (of ',length(teff.non),'non-member stars)',sep=""))
        legend("topleft",legend=leg.text,inset=0.02,pch=c(1,3),col=c('blue','red'))
      } else {
        num.stars <- sum(!is.na(teff.members))
        leg.text <- c(paste(num.stars,' (of ',length(teff.members),' member stars)',sep=""))
        legend("topleft",legend=leg.text,inset=0.02,pch=1,col='blue')
      }
      if (sum(!is.na(teff.members)) != 0) {
        mean.feh.members <- mean(feh.members,na.rm=TRUE)
        sd.feh.members <- sd(feh.members,na.rm=TRUE)
        title(main=bquote(paste('[Fe/H] = (', .(round(mean.feh.members, digits=2)) %+-% .(round(sd.feh.members,digits=2)),')',sep="" )))
        #            
        lin.reg <- lm(feh.members ~ teff.members)
        x.line <- seq(1000,10000,100)
        y.corr <- lin.reg$coefficients[1]+lin.reg$coefficients[2]*x.line
        lines(x.line,y.corr,lwd=2,lty=1,col='red')
      }
      if (sum(!is.na(teff.members)) >= 3) {
        #                legend("bottom", legend=bquote(paste('Y = (', .(round(lin.reg$coefficients[1], digits=2)) %+-% .(round(coef(summary(lin.reg))[1,"Std. Error"],digits=2)), ') + (',  .(round(lin.reg$coefficients[2], digits=2)) %+-% .(round(coef(summary(lin.reg))[2,"Std. Error"],digits=2)), ')*X',sep="" )), lty=1, lwd=2, col='red', cex=0.8)
      } else {
        #                legend("bottom", legend=bquote(paste('Y = (', .(round(lin.reg$coefficients[1], digits=2)), ') + (',  .(round(lin.reg$coefficients[2], digits=2)), ')*X',sep="" )), lty=1, lwd=2, col='red', cex=0.8)
      }
      #
      axis(side=1, at=seq(2000,10000,500), labels=seq(2000,10000,500), cex.axis=1.5, lwd=2)
      axis(side=1, at=seq(2000,10000,100), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
      axis(side=2, at=seq(-5.0,10.0,0.5), labels=seq(-5.0,10.0,0.5), cex.axis=1.5, lwd=2)
      axis(side=2, at=seq(-5.0,10.0,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
      #
      plotCI(logg.members,feh.members, ylim=range.feh, xlim=range.logg, ylab="[Fe/H] (dex)", xlab=paste("log g (dex)"), pch=1, cex=2.5, lwd=3, cex.lab=1.7, cex.axis=1.7, axes=FALSE, cex.main=1.5, xaxt='n', yaxt='n', col="blue")
      box(lwd=2)
      if (length(teff.non) != 0) {
        points(logg.non,feh.non, pch=3, col="red", cex=1.5, lwd=2)
        num.stars <- sum(!is.na(teff.members))
        num.non <- sum(!is.na(teff.non))
        leg.text <- c(paste(num.stars,' (of ',length(teff.members),' member stars)',sep=""),paste(num.non,' (of ',length(teff.non),'non-member stars)',sep=""))
        legend("topleft",legend=leg.text,inset=0.02,pch=c(1,3),col=c('blue','red'))
      } else {
        num.stars <- sum(!is.na(teff.members))
        leg.text <- c(paste(num.stars,' (of ',length(teff.members),' member stars)',sep=""))
        legend("topleft",legend=leg.text,inset=0.02,pch=1,col='blue')
      }
      if (sum(!is.na(teff.members)) != 0) {
        mean.feh.members <- mean(feh.members,na.rm=TRUE)
        sd.feh.members <- sd(feh.members,na.rm=TRUE)
        title(main=bquote(paste('[Fe/H] = (', .(round(mean.feh.members, digits=2)) %+-% .(round(sd.feh.members,digits=2)),')',sep="" )))
        #
        lin.reg <- lm(feh.members ~ logg.members)
        x.line <- seq(-2.0,8.0,0.5)
        y.corr <- lin.reg$coefficients[1]+lin.reg$coefficients[2]*x.line
        lines(x.line,y.corr,lwd=2,lty=1,col='red')
      }
      if (sum(!is.na(teff.members)) >= 3) {
        #                legend("bottom", legend=bquote(paste('Y = (', .(round(lin.reg$coefficients[1], digits=2)) %+-% .(round(coef(summary(lin.reg))[1,"Std. Error"],digits=2)), ') + (',  .(round(lin.reg$coefficients[2], digits=2)) %+-% .(round(coef(summary(lin.reg))[2,"Std. Error"],digits=2)), ')*X',sep="" )), lty=1, lwd=2, col='red',cex=0.8)
      } else {
        #                legend("bottom", legend=bquote(paste('Y = (', .(round(lin.reg$coefficients[1], digits=2)), ') + (',  .(round(lin.reg$coefficients[2], digits=2)), ')*X',sep="" )), lty=1, lwd=2, col='red',cex=0.8)
      }
      #
      axis(side=1, at=seq(-1.0,6.0,0.5), labels=seq(-1.0,6.0,0.5), cex.axis=1.5, lwd=2)
      axis(side=1, at=seq(-1.0,6.0,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
      axis(side=2, at=seq(-5.0,10.0,0.5), labels=seq(-5.0,10.0,0.5), cex.axis=1.5, lwd=2)
      axis(side=2, at=seq(-5.0,10.0,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
      dev.off()
    }
  }
  #
  for (each.cluster in names.clusters) {
    for (each.setup in c("U520","U580")) {
      filter.cluster <- ((in.nodes.ids[,col.gesfld] == each.cluster) & (grepl(each.setup,in.nodes.ids[,col.setup])))
      if (sum(filter.cluster) != 0) { # There are stars in that cluster observed with that setup
        if (sum(is.na(in.nodes.ids[filter.cluster,col.vel])) != length(in.nodes.ids[filter.cluster,col.vel])) { # There are some RVs there
          mean.rv <- mean(in.nodes.ids[filter.cluster,col.vel],na.rm=TRUE)
          sd.rv <- sd(in.nodes.ids[filter.cluster,col.vel],na.rm=TRUE)
        } else { # Else, everybody is a member
          mean.rv = 0.0
          sd.rv = 10000.0
        }
        if (length(dim(in.nodes.param)) == 3) {
          for (each.node in dimnames(in.nodes.param)[[2]]) {
            node.col <- which(dimnames(in.nodes.param)[[2]] == each.node)
            plot.the.cluster(each.cluster)
          }
        } else {
          for (each.node in c('Recommended_WG11')) {
            plot.the.cluster(each.cluster)
          }
        }
      }
    }
  }
}
###############################################
load.ges.dr5.data <- function(columns.abun,columns.metadata,path.to.files='/work/rsmiljanic/Survey/DR5/Support/',abun.file='gesiDR5RecommendedAstroAnalysis.fits',metadata.file='results13_11_57_46_531.fits') {
  # Read file with the abundances
  abun.idr5 <- readFITS(file=paste0(oath.to.files,abun.file), hdu=1, maxLines=25000)
  num.stars <- length(abun.idr5$col[[2]])
  # Extract the data in the columns listed in columns.abun
  extr.abun.idr5 <- matrix(' ',num.stars,length(columns.abun))
  colnames(extr.abun.idr5) <- columns.abun
  extr.abun.idr5 <- as.data.frame(extr.abun.idr5)
  #
  for (ik in seq(1,length(columns.abun),1)) {
    numb.col <- which(abun.idr5$colNames == columns.abun[ik])
    extr.abun.idr5[,ik] <- str_trim(abun.idr5$col[[numb.col]])
  }
  ges.dr5 <- extr.abun.idr5
  # Convert some columns to numeric
  numeric.columns <- columns.abun[!(columns.abun %in% c("CNAME", "GESFIELD", "GESOBJECT", "CONSTFILES", "GESTYPE", "RECGRATINGS", "RECWG"))]
  #
  for (each.col in numeric.columns) {
    this.col <- which(colnames(ges.dr5) == each.col)
    ges.dr5[,this.col] <- as.numeric(as.vector(ges.dr5[,this.col]))
    ges.dr5[(ges.dr5[,this.col] < -9000),this.col] <- NA
  }
  
  # Read file with metadata
  dr5.target <- readFITS(file=paste0(path.to.files,metadata.file), hdu=1, maxLines=55000)
  num.stars <- length(dr5.target$col[[2]])
  extr.dr5.target <- matrix(' ',num.stars,length(columns.metadata))
  colnames(extr.dr5.target) <- columns.metadata
  extr.dr5.target <- as.data.frame(extr.dr5.target)
  
  for (ik in seq(1,length(columns.metadata),1)) {
    numb.col <- which(dr5.target$colNames == columns.metadata[ik])
    extr.dr5.target[,ik] <- str_trim(dr5.target$col[[numb.col]])
  }
  # Convert some columns to numeric
  numeric.columns <- columns.metadata[!(columns.metadata %in% c("CNAME", "GESFIELD", "GESOBJECT", "CONSTFILES", "GESTYPE", "RECGRATINGS", "RECWG"))]
  #
  for (each.col in numeric.columns) {
    this.col <- which(colnames(extr.dr5.target) == each.col)
    extr.dr5.target[,this.col] <- as.numeric(as.vector(extr.dr5.target[,this.col]))
    extr.dr5.target[(extr.dr5.target[,this.col] < -9000),this.col] <- NA
  }
  
  
  #Join both tables
  order.abun <- order(ges.dr5$CNAME)
  order.target <- order(extr.dr5.target$CNAME)
  ges.dr5 <- ges.dr5[order.abun,]
  extr.dr5.target <- extr.dr5.target[order.target,]
  from.ges <- (ges.dr5$CNAME %in% extr.dr5.target$CNAME)
  from.target <- (extr.dr5.target$CNAME %in% ges.dr5$CNAME)
  ges.dr5 <- ges.dr5[from.ges,]
  extr.dr5.target <- extr.dr5.target[from.target,]
  drops <- c('CNAME')
  extr.dr5.target <- extr.dr5.target[,(!names(extr.dr5.target) %in% drops)]
  ges.dr5 <- cbind(ges.dr5,extr.dr5.target)
  
  #Return ges.dr5
  return(ges.dr5)
}

node.trends.with.recom <- function(in.nodes.param,recom.param,in.nodes.flags,release=substr(this.release,2,4),dir.to.save.plots='/work/rsmiljanic/Survey/DR6/Codes/Plots/') {
  #
  mean.sd.deltas <- matrix(NA,nrow=length(list.nodes),ncol=6)
  recom.param$TEFF[(recom.param$TEFF > 8000)] <- NA
  for (each.node in list.nodes) {
    ik.node <- which(list.nodes == each.node)
    teff.plotname <- paste(dir.to.save.plots,'wg11_',release,'_delta_teff_',each.node,'_recom.eps',sep="")
    logg.plotname <- paste(dir.to.save.plots,'wg11_',release,'_delta_logg_',each.node,'_recom.eps',sep="")
    feh.plotname <-  paste(dir.to.save.plots,'wg11_',release,'_delta_feh_',each.node,'_recom.eps',sep="")
    teff.vector <- as.numeric(as.vector(in.nodes.param[,ik.node,1]))
    logg.vector <- as.numeric(as.vector(in.nodes.param[,ik.node,3]))
    feh.vector <- as.numeric(as.vector(in.nodes.param[,ik.node,5]))
    delta.teff <- (teff.vector-recom.param$TEFF)
    delta.logg <- (logg.vector-recom.param$LOGG)
    delta.feh <- (feh.vector-recom.param$FEH)
    mean.sd.deltas[ik.node,1] <- mean(delta.teff,na.rm=T)
    mean.sd.deltas[ik.node,2] <- sd(delta.teff,na.rm=T)
    mean.sd.deltas[ik.node,3] <- mean(delta.logg,na.rm=T)
    mean.sd.deltas[ik.node,4] <- sd(delta.logg,na.rm=T)
    mean.sd.deltas[ik.node,5] <- mean(delta.feh,na.rm=T)
    mean.sd.deltas[ik.node,6] <- sd(delta.feh,na.rm=T)
    #
    postscript(file=teff.plotname,horizontal=FALSE, onefile=FALSE, height=6, width=7, pointsize=10)
    par(oma=c(0.5,0.5,0.8,0.5))
    par(mar=c(5,6,2,2))
    par(mgp=c(4,1.2,0))
    plot.title <- paste('UVES FGK-type: ',release,' ',each.node, sep=" ")
    plot(recom.param$TEFF,delta.teff,xlim=c(3500,7000), ylab=bquote(Delta ~ 'Teff (Node-Recom)'), xlab=expression(paste("Recommended ",T[eff], " (K)")), pch="*", cex=2.5, cex.lab=1.7, cex.axis=1.7, axes=FALSE, xaxt='n', yaxt='n', col="black", main=plot.title)
    abline(h=0,col='red',lwd=3)
    box(which="plot", lwd=3)
    axis(side=1, at=seq(3000,7500,500), labels=seq(3000,7500,500), cex.axis=1.5, lwd=2)
    axis(side=1, at=seq(3000,7500,100), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    axis(side=2, at=seq(-5000,5000,500), labels=seq(-5000,5000,500), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(-5000,5000,100), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    mean.teff <- mean(delta.teff,na.rm = T)
    sd.teff <- sd(delta.teff,na.rm = T)
    abline(h=mean.teff,col='blue',lwd=3)
    abline(h=mean.teff+2*sd.teff,col='blue',lwd=3,lty=2)
    abline(h=mean.teff-2*sd.teff,col='blue',lwd=3,lty=2)
    legend("topright", legend=bquote(paste('Mean = (', .(round(mean.teff, digits=2)) %+-% .(round(sd.teff,digits=2)),' )  ',sep="" )), lty=1, lwd=2, col='blue', cex=0.7)
    #legend("topright", legend=bquote(paste(Delta, ' = (', .(round(lin.reg$coefficients[1], digits=2)) %+-% .(round(coef(summary(lin.reg))[1,"Std. Error"],digits=2)), ') + (',  .(round(lin.reg$coefficients[2], digits=2)) %+-% .(round(coef(summary(lin.reg))[2,"Std. Error"],digits=2)), ')*',T[eff],'(Ref)',sep="" )), lty=1, lwd=2, col='red', cex=0.7)
    dev.off()
    #
    postscript(file=logg.plotname,horizontal=FALSE, onefile=FALSE, height=6, width=7, pointsize=10)
    par(oma=c(0.5,0.5,0.8,0.5))
    par(mar=c(5,6,2,2))
    par(mgp=c(4,1.2,0))
    plot.title <- paste('UVES FGK-type: ',release,' ',each.node, sep=" ")
    plot(recom.param$LOGG,delta.logg,xlim=c(0.0,5.5), ylab=bquote(Delta ~ 'log g (Node-Recom)'), xlab=expression(paste("Recommended ",log, "g (dex)")), pch="*", cex=2.5, cex.lab=1.7, cex.axis=1.7, axes=FALSE, xaxt='n', yaxt='n', col="black", main=plot.title)
    abline(h=0,col='red',lwd=3)
    box(which="plot", lwd=3)
    axis(side=1, at=seq(-7,7,0.5), labels=seq(-7,7,0.5), cex.axis=1.5, lwd=2)
    axis(side=1, at=seq(-7,7,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    axis(side=2, at=seq(-7,7,0.5), labels=seq(-7,7,0.5), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(-7,7,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    mean.logg <- mean(delta.logg,na.rm = T)
    sd.logg <- sd(delta.logg,na.rm = T)
    abline(h=mean.logg,col='blue',lwd=3)
    abline(h=mean.logg+2*sd.logg,col='blue',lwd=3,lty=2)
    abline(h=mean.logg-2*sd.logg,col='blue',lwd=3,lty=2)
    legend("topright", legend=bquote(paste('Mean = (', .(round(mean.logg, digits=2)) %+-% .(round(sd.logg,digits=2)),' )  ',sep="" )), lty=1, lwd=2, col='blue', cex=0.7)
    
    dev.off()      
    #
    postscript(file=feh.plotname,horizontal=FALSE, onefile=FALSE, height=6, width=7, pointsize=10)
    par(oma=c(0.5,0.5,0.8,0.5))
    par(mar=c(5,6,2,2))
    par(mgp=c(4,1.2,0))
    plot.title <- paste('UVES FGK-type: ',release,' ',each.node, sep=" ")
    plot(recom.param$FEH,delta.feh,xlim=c(-3.5,1.0), ylab=bquote(Delta ~ '[Fe/H] (Node-Recom)'), xlab=expression(paste("Recommended ","[Fe/H] (dex)")), pch="*", cex=2.5, cex.lab=1.7, cex.axis=1.7, axes=FALSE, xaxt='n', yaxt='n', col="black", main=plot.title)
    abline(h=0,col='red',lwd=3)
    box(which="plot", lwd=3)
    axis(side=1, at=seq(-5,5,0.5), labels=seq(-5,5,0.5), cex.axis=1.5, lwd=2)
    axis(side=1, at=seq(-5,5,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    axis(side=2, at=seq(-5,5,0.5), labels=seq(-5,5,0.5), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(-5,5,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    mean.feh <- mean(delta.feh,na.rm = T)
    sd.feh <- sd(delta.feh,na.rm = T)
    abline(h=mean.feh,col='blue',lwd=3)
    abline(h=mean.feh+2*sd.feh,col='blue',lwd=3,lty=2)
    abline(h=mean.feh-2*sd.feh,col='blue',lwd=3,lty=2)
    legend("topright", legend=bquote(paste('Mean = (', .(round(mean.feh, digits=2)) %+-% .(round(sd.feh,digits=2)),' )  ',sep="" )), lty=1, lwd=2, col='blue', cex=0.7)
    
    dev.off()      
    #
    postscript(file=paste(dir.to.save.plots,'wg11_',release,'_where_',each.node,'_failed_teff.eps',sep=""),horizontal=FALSE, onefile=FALSE, height=6, width=7, pointsize=10)
    par(oma=c(0.5,0.5,0.8,0.5))
    par(mar=c(5,6,2,2))
    par(mgp=c(4,1.2,0))
    plot.title <- paste('UVES FGK-type: ',release,' ',each.node, sep=" ")
    filter.failed.teff <- !is.na(recom.param$TEFF) & is.na(teff.vector) & (recom.param$TEFF >= 3000) & (recom.param$TEFF <= 7000) & !grepl('10399',nodes.flags[,ik.node,3])
    test <- hist(recom.param$TEFF[filter.failed.teff],breaks=seq(3000,7500,500),xlim=c(3000,7000), ylab='Number', xlab='Teff (K)', cex=2.5, cex.lab=1.7, cex.axis=1.7, axes=FALSE, xaxt='n', yaxt='n', main=plot.title)
    hist(recom.param$TEFF[filter.failed.teff],breaks=seq(3000,7500,500),xlim=c(3000,7000), ylab='Number', xlab='Teff (K)', cex=2.5, cex.lab=1.7, cex.axis=1.7, axes=FALSE, xaxt='n', yaxt='n', main=plot.title)
    text(4000,y=max(test$counts),labels=bquote(.(sum(filter.failed.teff)) ~ 'stars'),cex=1.5)
    box(which="plot", lwd=3)
    axis(side=1, at=seq(3000,8000,500), labels=seq(3000,8000,500), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(0,2000,20), labels=seq(0,2000,20), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(0,2000,5), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    dev.off()
    #
    postscript(file=paste(dir.to.save.plots,'wg11_',release,'_where_',each.node,'_failed_logg.eps',sep=""),horizontal=FALSE, onefile=FALSE, height=6, width=7, pointsize=10)
    par(oma=c(0.5,0.5,0.8,0.5))
    par(mar=c(5,6,2,2))
    par(mgp=c(4,1.2,0))
    plot.title <- paste('UVES FGK-type: ',release,' ',each.node, sep=" ")
    filter.failed.logg <- !is.na(recom.param$LOGG) & is.na(logg.vector) & (recom.param$LOGG >= 0.0) & (recom.param$LOGG <= 5.0) & !grepl('10399',nodes.flags[,ik.node,3])
    test <- hist(recom.param$LOGG[filter.failed.logg],breaks=seq(0.0,5.0,0.2),xlim=c(0.0,5.0), ylab='Number', xlab='log g (dex)', cex=2.5, cex.lab=1.7, cex.axis=1.7, axes=FALSE, xaxt='n', yaxt='n', main=plot.title)
    hist(recom.param$LOGG[filter.failed.logg],breaks=seq(0.0,5.0,0.2),xlim=c(0.0,5.0), ylab='Number', xlab='log g (dex)', cex=2.5, cex.lab=1.7, cex.axis=1.7, axes=FALSE, xaxt='n', yaxt='n', main=plot.title)
    text(1.0,y=max(test$counts),labels=bquote(.(sum(filter.failed.logg)) ~ 'stars'),cex=1.5)
    box(which="plot", lwd=3)
    axis(side=1, at=seq(0.0,5.0,0.2), labels=seq(0.0,5.0,0.2), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(0,2000,20), labels=seq(0,2000,20), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(0,2000,5), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    dev.off()
    #
    postscript(file=paste(dir.to.save.plots,'wg11_',release,'_where_',each.node,'_failed_feh.eps',sep=""),horizontal=FALSE, onefile=FALSE, height=6, width=7, pointsize=10)
    par(oma=c(0.5,0.5,0.8,0.5))
    par(mar=c(5,6,2,2))
    par(mgp=c(4,1.2,0))
    plot.title <- paste('UVES FGK-type: ',release,' ',each.node, sep=" ")
    filter.failed.feh <- !is.na(recom.param$FEH) & is.na(feh.vector) & (recom.param$FEH >= -3.0) & (recom.param$FEH <= 1.0) & !grepl('10399',nodes.flags[,ik.node,3])
    test <- hist(recom.param$FEH[filter.failed.feh],breaks=seq(-3.0,1.0,0.2),xlim=c(-3.0,1.0), ylab='Number', xlab='[Fe/H] (dex)', cex=2.5, cex.lab=1.7, cex.axis=1.7, axes=FALSE, xaxt='n', yaxt='n', main=plot.title)
    hist(recom.param$FEH[filter.failed.feh],breaks=seq(-3.0,1.0,0.2),xlim=c(-3.0,1.0), ylab='Number', xlab='[Fe/H] (dex)', cex=2.5, cex.lab=1.7, cex.axis=1.7, axes=FALSE, xaxt='n', yaxt='n', main=plot.title)
    text(-2.0,y=max(test$counts),labels=bquote(.(sum(filter.failed.feh)) ~ 'stars'),cex=1.5)
    box(which="plot", lwd=3)
    axis(side=1, at=seq(-3.0,1.0,0.2), labels=seq(-3.0,1.0,0.2), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(0,2000,20), labels=seq(0,2000,20), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(0,2000,5), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    dev.off()
    #
    postscript(file=paste(dir.to.save.plots,'wg11_',release,'_where_',each.node,'_failed_snr.eps',sep=""),horizontal=FALSE, onefile=FALSE, height=6, width=7, pointsize=10)
    par(oma=c(0.5,0.5,0.8,0.5))
    par(mar=c(5,6,2,2))
    par(mgp=c(4,1.2,0))
    snr.vector <- as.numeric(as.vector(recom.param$SNR))
    plot.title <- paste('UVES FGK-type: ',release,' ',each.node, sep=" ")
    filter.failed.feh <- !is.na(recom.param$FEH) & is.na(feh.vector) & (recom.param$FEH >= -3.0) & (recom.param$FEH <= 1.0) & !grepl('10399',nodes.flags[,ik.node,3]) & (recom.param$SNR <= 1100)
    test <- hist(recom.param$SNR[filter.failed.feh],breaks=seq(0,1100,20),xlim=c(0,500), ylab='Number', xlab='SNR', cex=2.5, cex.lab=1.7, cex.axis=1.7, axes=FALSE, xaxt='n', yaxt='n', main=plot.title)
    hist(recom.param$SNR[filter.failed.feh],breaks=seq(0,1100,20),xlim=c(0,500), ylab='Number', xlab='SNR', cex=2.5, cex.lab=1.7, cex.axis=1.7, axes=FALSE, xaxt='n', yaxt='n', main=plot.title)
    text(250,y=max(test$counts),labels=bquote(.(sum(filter.failed.feh)) ~ 'stars'),cex=1.5)
    box(which="plot", lwd=3)
    axis(side=1, at=seq(0,2000,20), labels=seq(0,2000,20), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(0,2000,20), labels=seq(0,2000,20), cex.axis=1.5, lwd=2)
    axis(side=2, at=seq(0,2000,5), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
    dev.off()
  }
  return(mean.sd.deltas)
}



