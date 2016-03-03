#***********************************************************************************************************
#***********************************************************************************************************
#********* blend3.R
#*********
#********* Copyright (C) 2014 Diamond Light Source & Imperial College London
#*********
#********* Authors: James Foadi & Gwyndaf Evans
#*********
#********* This code is distributed under the BSD license, a copy of which is
#********* included in the root directory of this package.
#************************************************************************************************************
#************************************************************************************************************
# Third R module of program BLEND (later to be replaced by C++ code)
# It relies on files produced by blend.cpp and blend1.R.
#
#


#
# Function to find out batches list in an mtz file
mtz_batch_list <- function(hklin,rwin=FALSE)
{
 # Build keywords file
 linea <- sprintf("END")
 cat(linea,file="mtzdump_keywords.dat")

 # Run mtzdmp and save output
 if (rwin)
 {
  exemtzdmp <- shell(paste("mtzdump.exe hklin ",hklin," < mtzdump_keywords.dat",sep=""),intern=TRUE)
 }
 else
 {
  exemtzdmp <- system(paste("mtzdump hklin ",hklin," < mtzdump_keywords.dat",sep=""),intern=TRUE)
 }

 # Extract list of batches from whole mtzdmp output
 ctmp <- grep("Batch number:",exemtzdmp,fixed=TRUE)+1
 blist <- c()
 for (i in ctmp)
 {
  stringa <- exemtzdmp[i]
  ltmp <- strsplit(stringa," ")
  idx <- which(nchar(ltmp[[1]]) != 0)
  blist <- c(blist,as.integer(ltmp[[1]][idx][1]))
 }

 # Delete keywords file
 if (file.exists("mtzdump_keywords.dat")) emptyc <- file.remove("mtzdump_keywords.dat")
 
 return(blist)
}

#
# Run pointless to find Laue group and merge files (after alternate indexing)
merge_mtzs <- function(mtz_list,selection,mtzout,pointless_keys,hklref,rwin=FALSE)
{
 # Cycle through each file in the list and do the following:
 #    1) extract batch information from each mtz file;
 #    2) add EXCLUDE BATCH for each file needed, according to selection
 #    3) include appropriate HKLREF for alternate indexing;
 #    4) glue all mtz together into a single mtz

 # Form POINTLESS keyword file
 linea <- sprintf("TITLE Prepare file for multicrystal merging\n")
 cat(linea,file="pointless_keywords.dat")
 linea <- sprintf("NAME PROJECT xxx CRYSTAL yyy DATASET zzz\n")
 cat(linea,file="pointless_keywords.dat",append=TRUE)
 for (imtz in 1:length(mtz_list))
 {
  hklin <- mtz_list[imtz]
  linea <- sprintf("HKLIN %s\n",hklin)
  cat(linea,file="pointless_keywords.dat",append=TRUE)

  # Extract batch information from each mtz file
  blist <- mtz_batch_list(hklin,rwin=rwin)

  # Images to remove (if requested)
  itmp <- selection[imtz]+1
  if (itmp > blist[length(blist)]) inibatch <- -1
  if (itmp <= blist[length(blist)]) inibatch <- itmp
  finbatch <- blist[length(blist)]
  if (inibatch > 0)
  {
   if (length(mtz_list) > 1) linea <- sprintf("EXCLUDE FILE %d BATCH %d TO %d\n",imtz,inibatch,finbatch)
   if (length(mtz_list) == 1) linea <- sprintf("EXCLUDE BATCH %d TO %d\n",inibatch,finbatch)
   cat(linea,file="pointless_keywords.dat",append=TRUE)
  }
 }
 linea <- sprintf("HKLREF %s\n",hklref)
 cat(linea,file="pointless_keywords.dat",append=TRUE)
 linea <- sprintf("HKLOUT %s\n",mtzout)
 cat(linea,file="pointless_keywords.dat",append=TRUE)

 # Additional keywords, introduced by user through BLEND_KEYWORDS.dat file
 for (line in pointless_keys) cat(paste(line,"\n",sep=""),file="pointless_keywords.dat",append=TRUE)

 # End of keywords file
 linea <- sprintf("END\n")
 cat(linea,file="pointless_keywords.dat",append=TRUE)
 
 # Run program
 if (rwin)
 {
  stringa <- "pointless.exe < pointless_keywords.dat"
  exepointless <- shell(stringa,intern=TRUE,ignore.stderr=TRUE)
 }
 else
 {
  stringa <- "pointless < pointless_keywords.dat"
  exepointless <- system(stringa,intern=TRUE,ignore.stderr=TRUE)
 }
 
 # Delete keywords file
 if (file.exists("pointless_keywords.dat")) emptyc <- file.remove("pointless_keywords.dat")

 return(exepointless)
}

#
# Run AIMLESS on specified datasets selection and collect statistical information in a data frame.
merge_datasets <- function(mtz_names,selection,suffix,pointless_keys,aimless_keys,resomin=NULL,resomax=NULL,nref=1,rwin=FALSE)
{
 # Computes merging statistics for all couples of datasets

 # Read content of file mtz_names
 if (file.exists(mtz_names)) indata <- read.table(file=mtz_names)
 if (!file.exists(mtz_names)) stop("Function argument mtz_names does not correspond to any valid existing file")

 # Turn factors into characters
 indata[,1] <- as.character(indata[,1])

 # Has CCP4 been set up?
 ccp4 <- Sys.getenv("CCP4")
 if (nchar(ccp4) == 0) stop("You need to set up environment for running ccp4 programs")


 # Reference file for indexing
 idxref <- as.integer(nref)
 if (is.na(idxref)) hklref <- nref
 if (!is.na(idxref)) hklref <- indata[idxref,1]
 

 # Run POINTLESS and AIMLESS
 mtz_list <- indata[selection,1]
 sele <- indata[selection,3]

 # Use copy of reference data set in case it belongs to group
 if (hklref %in% mtz_list)
 {
  file.copy(from=hklref,to="copy_of_ref_file.mtz",overwrite=TRUE)
  hklref <- "copy_of_ref_file.mtz"
 }
 cat("Collating multiple mtz into a single mtz ...\n")
 linea_in <- paste(suffix[1],paste("unscaled_",suffix[2],".mtz",sep=""),sep="/")
 exemerge <- merge_mtzs(mtz_list=mtz_list,selection=sele,mtzout=linea_in,pointless_keys=pointless_keys,hklref=hklref,rwin=rwin)
 if (file.exists("copy_of_ref_file.mtz")) file.remove("copy_of_ref_file.mtz")
 fatal_error <- grep("#-------------",exemerge,fixed=TRUE)
 if (length(fatal_error) == 0 &
     length(grep("FATAL ERROR message:",exemerge,fixed=TRUE)) == 0 &
     length(grep("Stopping",exemerge,fixed=TRUE)) == 0 &
     length(exemerge) != 0)
 {

  # Rename POINTLESS log
  log_file <- paste(suffix[1],paste("pointless_",suffix[2],".log",sep=""),sep="/")
  cat(exemerge,file=log_file,sep="\n")
  #for (linea in exemerge)
  #{
  # linea <- paste(linea,"\n",sep="")
  # cat(linea,file=log_file,append=TRUE)
  #}

  # Aimless keywords file
  linea_out <- paste(suffix[1],paste("scaled_",suffix[2],".mtz",sep=""),sep="/")
  linea <- "TITLE Scale multiple datasets\n"
  cat(linea,file="aimless_keywords.dat")
  linea <- sprintf("HKLIN %s\n",linea_in)
  cat(linea,file="aimless_keywords.dat",append=TRUE)
  linea <- sprintf("HKLOUT %s\n",linea_out)
  cat(linea,file="aimless_keywords.dat",append=TRUE)


  # Check RESOLUTION keyword is already assigned
  ncnd <- 0
  if (length(aimless_keys) > 0)         # i.e. if aimless_keys is different from NULL
  {
   tmp <- strsplit(aimless_keys," ")
   for (j in 1:length(tmp))
   {
    idx <- which(nchar(tmp[[j]]) != 0)
    cnd <- ((substr(tmp[[j]][idx],1,1) == "R") | (substr(tmp[[j]][idx],1,1) == "r")) &
           ((substr(tmp[[j]][idx],2,2) == "E") | (substr(tmp[[j]][idx],2,2) == "e")) &
           ((substr(tmp[[j]][idx],3,3) == "S") | (substr(tmp[[j]][idx],3,3) == "s")) &
           ((substr(tmp[[j]][idx],4,4) == "O") | (substr(tmp[[j]][idx],4,4) == "o"))
    ncnd <- ncnd+sum(cnd)
   }
  }
  if (length(aimless_keys) == 0) ncnd <- 0
  if (ncnd == 0)
  {
   if (!is.null(resomin) & !is.null(resomax))
   {
    linea <- sprintf("RESOLUTION LOW %f HIGH %f \n",resomin,resomax)
    cat(linea,file="aimless_keywords.dat",append=TRUE)
   }
   if (!is.null(resomin) & is.null(resomax))
   {
    linea <- sprintf("RESOLUTION LOW %f\n",resomin)
    cat(linea,file="aimless_keywords.dat",append=TRUE)
   }
   if (is.null(resomin) & !is.null(resomax))
   {
    linea <- sprintf("RESOLUTION HIGH %f\n",resomax)
    cat(linea,file="aimless_keywords.dat",append=TRUE)
   }
  }

  # Input by user through BLEND_KEYWORDS.dat
  for (line in aimless_keys) cat(paste(line,"\n",sep=""),file="aimless_keywords.dat",append=TRUE)

  # End of keywords file for AIMLESS
  linea <- sprintf("END\n")
  cat(linea,file="aimless_keywords.dat",append=TRUE)

  # Run AIMLESS
  cat("Running AIMLESS on the unscaled file ...\n")
  if (rwin)
  {
   stringa <- sprintf("aimless.exe < aimless_keywords.dat")
   exeaimless <- shell(stringa,intern=TRUE, ignore.stderr = TRUE)
  }
  else
  {
   stringa <- sprintf("aimless < aimless_keywords.dat")
   exeaimless <- system(stringa,intern=TRUE, ignore.stderr = TRUE)
  }

  # Collect and organize AIMLESS output
  nexeaimless <- length(exeaimless)
  if (length(grep("End of aimless job,", exeaimless[(nexeaimless - 1)], fixed = TRUE)) == 0)
  {
   Mstats <- data.frame(Rmeas=c(NA,NA,NA),Rpim=c(NA,NA,NA),Completeness=c(NA,NA,NA),Multiplicity=c(NA,NA,NA),
                        LowRes=c(NA,NA,NA),HighRes=c(NA,NA,NA))
   log_file <- paste(suffix[1],paste("aimless_",suffix[2],".log",sep=""),sep="/")
   for (linea in exeaimless)
   {
    linea <- paste(linea,"\n",sep="")
    cat(linea,file=log_file,append=TRUE)
   }
 
   # Clean directory from unnecessary files
   dircontents <- list.files("./")
   idx <- grep("ANOMPLOT",dircontents,fixed=TRUE)
   if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
   idx <- grep("CORRELPLOT",dircontents,fixed=TRUE)
   if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
   idx <- grep("NORMPLOT",dircontents,fixed=TRUE)
   if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
   idx <- grep("ROGUEPLOT",dircontents,fixed=TRUE)
   if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
   idx <- grep("ROGUES",dircontents,fixed=TRUE)
   if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
   idx <- grep("SCALES",dircontents,fixed=TRUE)
   if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
   idx <- grep("aimless_keywords",dircontents,fixed=TRUE)
   if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])

   # Remove files produced to merge mtz's
   dircontents <- list.files("./")
   idx <- grep("mtzdump",dircontents,fixed=TRUE)
   for (i in idx) file.remove(dircontents[i])
   idx <- grep("pointless",dircontents,fixed=TRUE)
   for (i in idx) file.remove(dircontents[i])
   #idx <- grep("reference",dircontents,fixed=TRUE)
   #for (i in idx) file.remove(dircontents[i])
 
   return(list(Mstats,exeaimless))
  }

  # Aimless output is OK. Carry on collecting information
  # Extract info on B factors and scales
  gRmeas <- grep("Rmeas (all",exeaimless,fixed=TRUE)
  lineRmeas <- exeaimless[gRmeas]
  tmpsplit <- strsplit(lineRmeas," ")
  idx <- which(nchar(tmpsplit[[1]]) != 0)
  Rmeas_values <- as.numeric(tmpsplit[[1]][idx[6:8]]) 
  gRpim <- grep("Rpim (all",exeaimless,fixed=TRUE)
  lineRpim <- exeaimless[gRpim]
  tmpsplit <- strsplit(lineRpim," ")
  idx <- which(nchar(tmpsplit[[1]]) != 0)
  Rpim_values <- as.numeric(tmpsplit[[1]][idx[6:8]]) 
  gCompleteness <- grep("Completeness  ",exeaimless,fixed=TRUE)
  lineCompleteness <- exeaimless[gCompleteness]
  tmpsplit <- strsplit(lineCompleteness," ")
  idx <- which(nchar(tmpsplit[[1]]) != 0)
  Completeness_values <- as.numeric(tmpsplit[[1]][idx[2:4]]) 
  gMultiplicity <- grep("Multiplicity  ",exeaimless,fixed=TRUE)
  lineMultiplicity <- exeaimless[gMultiplicity]
  tmpsplit <- strsplit(lineMultiplicity," ")
  idx <- which(nchar(tmpsplit[[1]]) != 0)
  Multiplicity_values <- as.numeric(tmpsplit[[1]][idx[2:4]]) 
  gLowResos <- grep("Low resolution limit",exeaimless,fixed=TRUE)
  lineLowResos <- exeaimless[gLowResos]
  tmpsplit <- strsplit(lineLowResos," ")
  idx <- which(nchar(tmpsplit[[1]]) != 0)
  LowResos_values <- as.numeric(tmpsplit[[1]][idx[4:6]])
  gHighResos <- grep("High resolution limit",exeaimless,fixed=TRUE)
  lineHighResos <- exeaimless[gHighResos]
  tmpsplit <- strsplit(lineHighResos," ")
  idx <- which(nchar(tmpsplit[[1]]) != 0)
  HighResos_values <- as.numeric(tmpsplit[[1]][idx[4:6]])
 
  Mstats <- data.frame(Rmeas=Rmeas_values,Rpim=Rpim_values,Completeness=Completeness_values,Multiplicity=Multiplicity_values,
                       LowRes=LowResos_values,HighRes=HighResos_values)

  # Write AIMLESS log file
  log_file <- paste(suffix[1],paste("aimless_",suffix[2],".log",sep=""),sep="/")
  cat(exeaimless,file=log_file,sep="\n")
  #for (linea in exeaimless)
  #{
  # linea <- paste(linea,"\n",sep="")
  # cat(linea,file=log_file,append=TRUE)
  #}

  # Remove files produced to merge mtz's
  dircontents <- list.files("./")
  idx <- grep("mtzdump",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (i in idx) file.remove(dircontents[i])
  idx <- grep("pointless",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (i in idx) file.remove(dircontents[i])
  #idx <- grep("reference",dircontents,fixed=TRUE)
  #if (length(idx) != 0) for (i in idx) file.remove(dircontents[i])
 
  # Remove files produced to scale mtz's
  dircontents <- list.files("./")
  idx <- grep("ANOMPLOT",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
  idx <- grep("CORRELPLOT",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
  idx <- grep("NORMPLOT",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
  idx <- grep("ROGUEPLOT",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
  idx <- grep("ROGUES",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
  idx <- grep("SCALES",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
  idx <- grep("aimless_keywords",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])

  return(list(Mstats,exeaimless))
 }
 else
 {
  Mstats <- data.frame(Rmeas=c(NA,NA,NA),Rpim=c(NA,NA,NA),Completeness=c(NA,NA,NA),Multiplicity=c(NA,NA,NA),
                       LowRes=c(NA,NA,NA),HighRes=c(NA,NA,NA))

  # Rename POINTLESS log
  log_file <- paste(suffix[1],paste("pointless_",suffix[2],".log",sep=""),sep="/")
  cat(exemerge,file=log_file,sep="\n")
  #for (linea in exemerge)
  #{
  # linea <- paste(linea,"\n",sep="")
  # cat(linea,file=log_file,append=TRUE)
  #}
  
  return(list(Mstats,exemerge))
 }
}

#
# Function to interpret input for the combination mode.
# 1) individual numbers, like 2 4 8 17 etc., simply means that the group is formed of individual data sets: "2 4 8 17"
# 2) individual numbers separated by commas, like 2,4,8,17, means: "2 4 8 17"
# 3) numbers separated by "-" means a range, like 12-16, means: "12 13 14 15 16"
# 4) combination of the above, like 1 3 4-9 12,16, means: "1 3 4 5 6 7 8 9 12 16"
# 5) one number within square bracket means add elements of a cluster. For example cluster 8 might include 
#    data sets 2, 7, 13, 17 and 23. Thus [8] means: "2 7 13 17 23"
# 6) one number within curly brackets means to subtract that number from a list. For example, using the same
#    cluster as before, the expression [8] {13,23} means "2 7 17"
# 7) same points 1 to 4 apply to everything inside a curly bracket.
interpretCombine <- function(args,clusters)
{
 # Number of character strings
 n <- length(args)

 # Loop over all character strings and extract information as individual numbers
 addList <- c()
 subtractList <- c()
 for (i in 1:n)
 {
  stmp <- args[i]

  # Differentiate between "no brackets" and brackets ("square" or "double square")
  m <- nchar(stmp)

  # Square brackets
  if (substr(stmp,1,1) == "[" & substr(stmp,nchar(stmp),nchar(stmp)) == "]" &
      substr(stmp,2,2) != "[" & substr(stmp,nchar(stmp)-1,nchar(stmp)-1) != "]")
  {
   gcore <- substr(stmp,2,nchar(stmp)-1)
   glist <- interpretGcore(gcore)

   # Expand list of cluster as clusters elements
   for (i in glist)
   {
    addList <- c(addList,clusters[[i]])
   }
  }

  # Double square brackets
  if (substr(stmp,1,1) == "[" & substr(stmp,2,2) == "[" &
      substr(stmp,nchar(stmp)-1,nchar(stmp)-1) == "]" & substr(stmp,nchar(stmp),nchar(stmp)) == "]")
  {
   gcore <- substr(stmp,3,nchar(stmp)-2)
   glist <- interpretGcore(gcore)
   subtractList <- c(subtractList,glist)
  }

  # No brackets
  if (substr(stmp,1,1) != "[" & substr(stmp,nchar(stmp),nchar(stmp)) != "]")
  {
   gcore <- stmp
   glist <- interpretGcore(gcore)
   addList <- c(addList,glist)
  }
 }
 addList <- unique(addList)
 if (length(subtractList) != 0) subtractList <- unique(subtractList)

 # Remove elements in subtractList from elements in addList
 if (length(subtractList) != 0) 
 {
  idx <- na.omit(match(subtractList,addList))
  if (length(idx) > 0) args <- addList[-idx]
  if (length(idx) == 0) args <- addList
 }
 if (length(subtractList) == 0) 
 {
  args <- addList
 }

 # Sorting before returning
 args <- sort(args)

 return(args)
}

interpretGcore <- function(gcore)
{
 # Find out number of blocks
 glist <- c()
 ltmp <- strsplit(gcore,",")
 if (length(ltmp[[1]]) > 1)
 {
  for (j in 1:length(ltmp[[1]]))
  {
   ltmp2 <- strsplit(ltmp[[1]][j],"-")
   if (length(ltmp2[[1]]) == 2)
   {
    if (as.integer(ltmp2[[1]][1]) <= as.integer(ltmp2[[1]][2])) glist <- c(glist,ltmp2[[1]][1]:ltmp2[[1]][2])
   }
   if (length(ltmp2[[1]]) == 1)
   {
    glist <- c(glist,as.integer(ltmp2[[1]][1]))
   }
  }
 }
 if (length(ltmp[[1]]) == 1)
 {
  ltmp2 <- strsplit(ltmp[[1]][1],"-")
  if (length(ltmp2[[1]]) == 2)
  {
   if (as.integer(ltmp2[[1]][1]) <= as.integer(ltmp2[[1]][2])) glist <- c(glist,ltmp2[[1]][1]:ltmp2[[1]][2])
  }
  if (length(ltmp2[[1]]) == 1)
  {
   glist <- c(glist,as.integer(ltmp2[[1]][1]))
  }
 }

 return(glist)
}

#
# Given a table "tbl", which is, in fact, a vector of strings obtained while reading
# a table from an ascii file, this function returns a list whose elements are vectors of
# values corresponding to the table columns, whose positions are given by vector "pns".
get_vectors_from_table <- function(tbl,pns)
{
 # Turn vector of strings into a list of vectors
 tmpsplit <- strsplit(tbl,split=" ",fixed=TRUE)
 m <- length(tmpsplit)
 ltmp <- list()
 for (i in 1:m)
 {
  vec <- tmpsplit[[i]]
  idx <- which(vec == "")
  ltmp <- c(ltmp,list(vec[-idx]))
 }
 tmpsplit <- lapply(ltmp,as.numeric)
 rm(ltmp)

 # Create final matrix
 n <- length(pns)
 M <- matrix(nrow=m,ncol=n)
 for (i in 1:m)
 {
  M[i,] <- tmpsplit[[i]][pns]
 }
 rm(tmpsplit)

 return(M)
}

#
# Pruning plan
pruning_plan <- function(tbl,completeness,cutf,icyc)
{
 # Partition all images in runs
 prt <- floor(tbl[,1]/1000)
 rns <- unique(prt)
 istart <- c()
 iend <- c()
 for (i in rns)
 {
  idx <- which(prt == i)
  istart <- c(istart,idx[1])
  iend <- c(iend,tail(idx,n=1))
 }
 Mruns <- matrix(c(istart,iend),ncol=2)

 # Number of images that need to be excluded (based on default or user value of "completeness")
 idx <- which(tbl[,3] >= completeness)
 n_im_elim <- length(idx)

 # Partition number of discarded images among all runs (proportionally to each dataset)
 neq <- c()
 for (irun in rns)
 {
  idx <- which(prt == irun)
  r <- n_im_elim*length(idx)/length(prt)
  neq <- c(neq,floor(r))
 }

 # Select run with highest mean Rmerge
 mean_Rmerge <- c()
 Images <- c()
 for (i in rns)
 {
  j <- i+1
  mean_Rmerge <- c(mean_Rmerge,mean(tbl[Mruns[j,1]:Mruns[j,2],2],na.rm=TRUE))
  Images <- c(Images,length(Mruns[j,1]:Mruns[j,2]))
 }
 mean_Rmerge[mean_Rmerge < 0] <- NA

 isel <- which(mean_Rmerge == max(mean_Rmerge,na.rm=TRUE))
 if (length(isel) == 0) isel <- 1

 # Find exact values of images to eliminate
 neq[isel] <- floor(cutf*neq[isel])
 if (neq[isel] > 0) ini_ima <- tbl[Mruns[isel,2]-neq[isel]+1,1]
 fin_ima <- tbl[Mruns[isel,2],1]
 if (neq[isel] == 0) ini_ima <- fin_ima
 
 # If number of images left after elimination is less than a hard value (5) suspend pruning
 HARD_VALUE <- 5
 idx <- tbl[Mruns[isel,1]:(Mruns[isel,2]-neq[isel]),1]
 if (length(idx) < HARD_VALUE)
 {
  neq[isel] <- 0
  stmp <- sprintf("Only %d images left for run %d.\n",length(idx),isel)
  cat(stmp)
 }

 return(list(ini_ima=ini_ima,fin_ima=fin_ima,n_im_elim=n_im_elim,neq=neq,isel=isel,mean_Rmerge=mean_Rmerge,Images=Images))
}

#
# Run AIMLESS on previously created unscaled file, with additional keywords for exclusion of images
simple_merge_datasets <- function(suffix,aimless_keys,rwin=FALSE)
{
 # AIMLESS keywords file
 cat(aimless_keys,file="aimless_keywords.dat",sep="\n")
 
 # Aimless scaled file
 linea_out <- paste(suffix[1],paste("scaled_",suffix[2],".mtz",sep=""),sep="/")

 # Run AIMLESS
 cat("Re-running AIMLESS on the unscaled file ...\n")
 if (rwin)
 {
  stringa <- sprintf("aimless.exe < aimless_keywords.dat")
  exeaimless <- shell(stringa,intern=TRUE, ignore.stderr = TRUE)
 }
 else
 {
  stringa <- sprintf("aimless < aimless_keywords.dat")
  exeaimless <- system(stringa,intern=TRUE, ignore.stderr = TRUE)
 }

 # Collect and organize AIMLESS output
 nexeaimless <- length(exeaimless)
 if (length(grep("End of aimless job,", exeaimless[(nexeaimless - 1)], fixed = TRUE)) == 0)
 {
  Mstats <- data.frame(Rmeas=c(NA,NA,NA),Rpim=c(NA,NA,NA),Completeness=c(NA,NA,NA),Multiplicity=c(NA,NA,NA),
                       LowRes=c(NA,NA,NA),HighRes=c(NA,NA,NA))
  log_file <- paste(suffix[1],paste("aimless_",suffix[2],".log",sep=""),sep="/")
  cat(exeaimless,file=log_file,sep="\n")
 
  # Clean directory from unnecessary files
  dircontents <- list.files("./")
  idx <- grep("ANOMPLOT",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
  idx <- grep("CORRELPLOT",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
  idx <- grep("NORMPLOT",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
  idx <- grep("ROGUEPLOT",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
  idx <- grep("ROGUES",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
  idx <- grep("SCALES",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
  idx <- grep("aimless_keywords",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])

  # Remove files produced to merge mtz's
  dircontents <- list.files("./")
  idx <- grep("mtzdump",dircontents,fixed=TRUE)
  for (i in idx) file.remove(dircontents[i])
 
  return(list(Mstats,exeaimless))
 }

 # Aimless output is OK. Carry on collecting information
 # Extract info on B factors and scales
 gRmeas <- grep("Rmeas (all",exeaimless,fixed=TRUE)
 lineRmeas <- exeaimless[gRmeas]
 tmpsplit <- strsplit(lineRmeas," ")
 idx <- which(nchar(tmpsplit[[1]]) != 0)
 Rmeas_values <- as.numeric(tmpsplit[[1]][idx[6:8]]) 
 gRpim <- grep("Rpim (all",exeaimless,fixed=TRUE)
 lineRpim <- exeaimless[gRpim]
 tmpsplit <- strsplit(lineRpim," ")
 idx <- which(nchar(tmpsplit[[1]]) != 0)
 Rpim_values <- as.numeric(tmpsplit[[1]][idx[6:8]]) 
 gCompleteness <- grep("Completeness  ",exeaimless,fixed=TRUE)
 lineCompleteness <- exeaimless[gCompleteness]
 tmpsplit <- strsplit(lineCompleteness," ")
 idx <- which(nchar(tmpsplit[[1]]) != 0)
 Completeness_values <- as.numeric(tmpsplit[[1]][idx[2:4]]) 
 gMultiplicity <- grep("Multiplicity  ",exeaimless,fixed=TRUE)
 lineMultiplicity <- exeaimless[gMultiplicity]
 tmpsplit <- strsplit(lineMultiplicity," ")
 idx <- which(nchar(tmpsplit[[1]]) != 0)
 Multiplicity_values <- as.numeric(tmpsplit[[1]][idx[2:4]]) 
 gLowResos <- grep("Low resolution limit",exeaimless,fixed=TRUE)
 lineLowResos <- exeaimless[gLowResos]
 tmpsplit <- strsplit(lineLowResos," ")
 idx <- which(nchar(tmpsplit[[1]]) != 0)
 LowResos_values <- as.numeric(tmpsplit[[1]][idx[4:6]])
 gHighResos <- grep("High resolution limit",exeaimless,fixed=TRUE)
 lineHighResos <- exeaimless[gHighResos]
 tmpsplit <- strsplit(lineHighResos," ")
 idx <- which(nchar(tmpsplit[[1]]) != 0)
 HighResos_values <- as.numeric(tmpsplit[[1]][idx[4:6]])
 
 Mstats <- data.frame(Rmeas=Rmeas_values,Rpim=Rpim_values,Completeness=Completeness_values,Multiplicity=Multiplicity_values,
                      LowRes=LowResos_values,HighRes=HighResos_values)

 # Write AIMLESS log file
 log_file <- paste(suffix[1],paste("aimless_",suffix[2],".log",sep=""),sep="/")
 cat(exeaimless,file=log_file,sep="\n")

 # Remove files produced to merge mtz's
 dircontents <- list.files("./")
 idx <- grep("mtzdump",dircontents,fixed=TRUE)
 if (length(idx) != 0) for (i in idx) file.remove(dircontents[i])
 
 # Remove files produced to scale mtz's
 dircontents <- list.files("./")
 idx <- grep("ANOMPLOT",dircontents,fixed=TRUE)
 if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
 idx <- grep("CORRELPLOT",dircontents,fixed=TRUE)
 if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
 idx <- grep("NORMPLOT",dircontents,fixed=TRUE)
 if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
 idx <- grep("ROGUEPLOT",dircontents,fixed=TRUE)
 if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
 idx <- grep("ROGUES",dircontents,fixed=TRUE)
 if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
 idx <- grep("SCALES",dircontents,fixed=TRUE)
 if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])
 idx <- grep("aimless_keywords",dircontents,fixed=TRUE)
 if (length(idx) != 0) for (jj in idx) file.remove(dircontents[jj])

 return(list(Mstats,exeaimless))
}


######################################################################################################################
######################################################################################################################

# Main

# Find out if running on Windows
ostuff <- Sys.info()
rwin = FALSE
if (ostuff["sysname"] == "Windows") rwin = TRUE

# Load content of previous R run
if (file.exists("BLEND.RData")) load("BLEND.RData",.GlobalEnv)
if (!file.exists("BLEND.RData"))
{
 msg <- "File BLEND.RData (produced by a previous run of BLEND in analysis mode) not found.\n"
 cat(msg)
 q(save = "no", status = 1, runLast = FALSE)
}

# To avoid warning messages set warn to a negative value 
options(warn = -1)

# Retrieve value from command line
args <- commandArgs(trailingOnly=TRUE)
combination_type <- as.integer(args[1])
tmp <- args
args <- tmp[2:length(tmp)]
rm(tmp)
cat("\n")
cat(paste("Input string:   ",paste(args,collapse=" ")))
n <- length(args)
cat("\n")
cat("Input request:\nData sets to be included:\n")
for (i in 1:n)
{
 stmp <- args[i]
 if (substr(stmp,1,1) != "[" & substr(stmp,nchar(stmp),nchar(stmp)) != "]")
 {
  gcore <- stmp
  glist <- interpretGcore(gcore)
  msg <- sprintf("                                              : %s\n",paste(glist,collapse=" "))
  cat(msg)
 }
}
cat("Input request:\nClusters with data sets to be included:\n")
for (i in 1:n)
{
 stmp <- args[i]
 if (substr(stmp,1,1) == "[" & substr(stmp,nchar(stmp),nchar(stmp)) == "]" &
     substr(stmp,2,2) != "[" & substr(stmp,nchar(stmp)-1,nchar(stmp)-1) != "]")
 {
  gcore <- substr(stmp,2,nchar(stmp)-1)
  glist <- interpretGcore(gcore)
  for (j in glist) 
  {
   if (j < 1 | j > length(groups[[1]]))
   {
    stmp <- sprintf(" USER ERROR! Cluster %3d does not exist.\n",j)
    cat(stmp)
    q(save = "no", status = 1, runLast = FALSE)
   }
   msg <- sprintf("         Cluster %3d - Includes data sets     : %s\n",j,paste(groups[[1]][[j]],collapse=" "))
   cat(msg)
  }
 }
}
cat("Input request:\nData sets to be excluded:\n")
for (i in 1:n)
{
 stmp <- args[i]
 if (substr(stmp,1,1) == "[" & substr(stmp,2,2) == "[" &
     substr(stmp,nchar(stmp)-1,nchar(stmp)-1) == "]" & substr(stmp,nchar(stmp),nchar(stmp)) == "]")
 {
  gcore <- substr(stmp,3,nchar(stmp)-2)
  glist <- interpretGcore(gcore)
  msg <- sprintf("                                              : %s\n",paste(glist,collapse=" "))
  cat(msg)
 }
}
cat("\n")
cat("Files to be combined:  ")
targs <- interpretCombine(args,groups[[1]])
if (length(targs) == 0) cat("None")
if (length(targs) > 0) cat(targs)
cat("\n")
cat("\n")
cln <- targs

# If user selection implies no data sets, return with warning message
if (length(cln) == 0)
{
 msg <- "WARNING! User selection of data sets for this combination implies no data sets.\n\n"
 cat(msg)
 q(save = "no", status = 0, runLast = FALSE)
}

# If one or more datasets have been discarded during the analysis mode, stop
tmp <- match(cln,maindf$cn)
if (sum(is.na(tmp)) > 0) stop("One or more input datasets have been discarded by BLEND during the analysis run, or have never been processed.")

# Find out whether the series of datasets in this exact order has been already processed
gidx <- c()                  # In case there's no GROUPS.info file
outdir <- "./combined_files"
if (file.exists(outdir))
{
 dircontents <- list.files(outdir)

 # Check whether number of files is OK in combined_files directory
 ltmp <- sapply(dircontents,FUN=strsplit,split=".",fixed=TRUE)
 tmp <- unlist(ltmp)
 lidx <- which(tmp == "log")
 midx <- which(tmp == "mtz")
 if (length(midx) != length(lidx))
 {
  messaggio <- "WARNING! Number of POINTLESS and AIMLESS logs in combined_files directory does not correspond to number of reflection files.\n\n"
  cat(messaggio)
 }

 # Extract serial number of groups already existing in file GROUPS.info
 mtmp <- paste(outdir,"/GROUPS.info",sep="")
 if (!file.exists(mtmp)) cat("WARNING! File GROUPS.info is missing in directory combined_files\n\n")
 if (file.exists(mtmp))
 {
  groups_info <- scan(mtmp,what="character",sep="\n",quiet=TRUE)
  gidx <- grep("Group",groups_info,fixed=TRUE)
  if (length(gidx) > 1)
  {
   ini <- 2
   fin <- c()
   for (i in 2:length(gidx))
   {
    fin <- c(fin,gidx[i]-1)
    ini <- c(ini,gidx[i]+1)
   }
   fin <- c(fin,length(groups_info))
  }
  if (length(gidx) == 1)
  {
   ini <- 2
   fin <- length(groups_info)
  }
  glist <- list()
  for (i in 1:length(ini))
  {
   gvec <- c()
   for (j in ini[i]:fin[i])
   {
    ltmp <- strsplit(groups_info[j]," ")
    idx <- which(nchar(ltmp[[1]]) != 0)
    gvec <- c(gvec,as.numeric(ltmp[[1]][idx[2]]))
   }
   glist <- c(glist,list(gvec))
  }

  # Check whether groups described in GROUPS.info matches files contained in combined_files
  if (length(lidx) != 2*length(glist)) cat("WARNING! Number of groups described in file GROUPS.info does not match number of files in directory combined_files.\n\n")
 }

 # Load content of MERGING_STATISTICS.info, if it exists
 mtmp <- paste(outdir,"/MERGING_STATISTICS.info",sep="")
 if (!file.exists(mtmp)) cat("WARNING! File MERGING_STATISTICS.info is missing in directory combined_files\n\n")
 if (file.exists(mtmp))
 {
  merging_statistics_info <- scan(mtmp,what="character",sep="\n",quiet=TRUE)
  idx <- grep("Rmeas",merging_statistics_info,fixed=TRUE)
  nn <- length(merging_statistics_info)-idx
  if (nn != length(gidx)) cat("WARNING! MERGING_STATISTICS.info has information contrasting with GROUPS_info.\n\n")
 }
}

# Find out resolutions corresponding to selected datasets
fdata <- read.table("FINAL_list_of_files.dat")
jdx <- match(cln,fdata$V2)
resosel <- max(fdata$V6[jdx])

# Extract keywords files for POINTLESS and AIMLESS
pointless_keys <- c()
aimless_keys <- c()
if (file.exists("BLEND_KEYWORDS.dat"))
{
 contents <- scan("BLEND_KEYWORDS.dat", what="character",sep="\n",quiet=TRUE)

 # BLEND KEYWORDS
 # Reference dataset (from BLEND KEYWORDS section)
 tmp <- grep("DATA",contents,fixed=TRUE)
 if (length(tmp) != 0)
 {
  tmp <- contents[tmp[1]]    # [1] in case somebody add DATAREF line more than once
  tmp <- strsplit(tmp," ")
  jdx <- which(nchar(tmp[[1]]) != 0)
  idxref_char <- tmp[[1]][jdx[length(jdx)]]
  idxref <- as.integer(idxref_char)
 } else idxref_char <- as.character(idxref)

 # Read in "FINAL_list_of_files.dat" to associate reference dataset to file path (just in some cases)
 indata <- read.table(file="FINAL_list_of_files.dat")
 indata[,1] <- as.character(indata[,1])
 tmpref <- as.integer(idxref_char)
 if (!is.na(tmpref)) idxref_char <- indata[idxref,1]
 if (!file.exists(idxref_char))
 {
  messaggio <- paste("Input reference dataset",idxref_char,"does not exist. Please input an existing file.\n")
  cat(messaggio)
  q(save = "no", status = 1, runLast = FALSE)
 }
 messaggio <- paste("Reference dataset used in case alternative indexing is needed: ",idxref_char,"\n",sep="")
 cat(messaggio)
 cat("\n")

 # Maximum number of cycles for BLEND -cP mode (from BLEND KEYWORDS section)
 tmp <- grep("MAXC",contents,fixed=TRUE)
 if (length(tmp) != 0)
 {
  tmp <- contents[tmp[1]]    # [1] in case somebody add MAXCYCLE line more than once
  tmp <- strsplit(tmp," ")
  jdx <- which(nchar(tmp[[1]]) != 0)
  maxcycle_char <- tmp[[1]][jdx[length(jdx)]]
  maxcycle <- as.integer(maxcycle_char)
 } else maxcycle <- 5

 # Minimum data completeness for BLEND -cP mode (from BLEND KEYWORDS section)
 tmp <- grep("COMP",contents,fixed=TRUE)
 if (length(tmp) != 0)
 {
  tmp <- contents[tmp[1]]    # [1] in case somebody add COMPLETENESS line more than once
  tmp <- strsplit(tmp," ")
  jdx <- which(nchar(tmp[[1]]) != 0)
  completeness_char <- tmp[[1]][jdx[length(jdx)]]
  completeness <- as.numeric(completeness_char)
 } else completeness <- 95.0

 # Fraction of images to be eliminated for BLEND -cp mode (from BLEND KEYWORDS section)
 tmp <- grep("CUTF",contents,fixed=TRUE)
 if (length(tmp) != 0)
 {
  tmp <- contents[tmp[1]]    # [1] in case somebody add COMPLETENESS line more than once
  tmp <- strsplit(tmp," ")
  jdx <- which(nchar(tmp[[1]]) != 0)
  cutfraction_char <- tmp[[1]][jdx[length(jdx)]]
  cutfraction <- as.numeric(cutfraction_char)
 } else cutfraction <- 0.5

 # POINTLESS KEYWORDS
 idxs <- grep("POINTLESS KEYWORDS",contents,fixed=TRUE)
 idxe <- grep("AIMLESS KEYWORDS",contents,fixed=TRUE)
 if ((idxe-idxs) > 1) for (line in contents[(idxs+1):(idxe-1)]) pointless_keys <- c(pointless_keys,line)

 # AIMLESS
 idxs <- idxe
 idxe <- length(contents)
 if ((idxe-idxs) > 0) for (line in contents[(idxs+1):idxe]) aimless_keys <- c(aimless_keys,line)
}

# Produce scaled dataset
if (!file.exists(outdir)) dir.create(outdir)   # Create "combined_files" directory, if it does not exist
ftail <- sprintf("%03d",(length(gidx)+1))
suffix <- c(outdir,ftail)
if (combination_type == 0)
{
 tmp <- merge_datasets("FINAL_list_of_files.dat",selection=cln,suffix,pointless_keys,aimless_keys,
                       resomax=resosel,nref=idxref_char,rwin=rwin)

 # Get table of average Rmerges per run
 if (!is.na(tmp[[1]]$Rmeas[1]))
 {
  # Read AIMLESS log and store some info
  log_aimless <- paste(suffix[1],paste("aimless_",suffix[2],".log",sep=""),sep="/")
  contents_aimless <- scan(log_aimless,what="character",sep="\n",quiet=TRUE)
 
  # Rmerge values
  gRmerge <- grep("N   Batch    Mn(I)   RMSdev  I/rms  Rmerge",contents_aimless,fixed=TRUE)
  istart <- gRmerge[1]+1
  iend <- gRmerge[2]-2
  line_Rmerge <- contents_aimless[istart:iend]
  #tbl <- get_vectors_from_table(line_Rmerge,c(2,6,9,13))
  tbl <- get_vectors_from_table(line_Rmerge,c(2,6,9))
  lista <- pruning_plan(tbl,100,1,0)
  mean_Rmerge <- lista$mean_Rmerge
  Images <- lista$Images
  rm(lista)

  # Output table for average Rmerge and Images per datasets plots 
  nruns_final <- length(mean_Rmerge)
  cat("\n")
  cat("$TABLE: Image number and average Rmerge for all datasets :\n")
  cat("$GRAPHS:    Number of images per dataset       : N : 1, 2 :\n")
  cat("       :    Average Rmerge per dataset         : N : 1, 3 :\n")
  cat("$$\n")
  cat("  Datasets   Images    mean(Rmerge)  $$\n")
  cat("$$\n")
  for (irun in 1:nruns_final)
  {
   stmp <- sprintf("%10d %8d %15.3f\n",cln[irun],Images[irun],mean_Rmerge[irun])
   cat(stmp)
  }
  cat("$$\n")
  cat("\n")
 }
}
if (combination_type == 1)
{
 cat("\n")
 cat("Working out strategy to improve merging statistics ...\n")
 cat("\n")

 # Create two lists to store aimless keywords and merging statistics.
 stored_results <- list()
 stored_aimkeys <- list()

 # Run default POINTLESS and AIMLESS to load starting information
 cat("### Cycle 0 ###\n")
 tmp <- merge_datasets("FINAL_list_of_files.dat",selection=cln,suffix,pointless_keys,aimless_keys,
                       resomax=resosel,nref=idxref_char,rwin=rwin)
 if (!is.na(tmp[[1]]$Rmeas[1]))
 {
  stored_results <- c(stored_results,list(tmp[[1]]))
  stored_aimkeys <- c(stored_aimkeys,list(aimless_keys))

  # Read AIMLESS log and store some info
  log_aimless <- paste(suffix[1],paste("aimless_",suffix[2],".log",sep=""),sep="/")
  contents_aimless <- scan(log_aimless,what="character",sep="\n",quiet=TRUE)

  # Rmeas, Rpim and Completeness for this (starting) cycle
  gRmeas <- grep("Rmeas (all",contents_aimless,fixed=TRUE)
  lineRmeas <- contents_aimless[gRmeas]
  tmpsplit <- strsplit(lineRmeas," ")
  idx <- which(nchar(tmpsplit[[1]]) != 0)
  Rmeas_values <- as.numeric(tmpsplit[[1]][idx[6:8]]) 
  gRpim <- grep("Rpim (all",contents_aimless,fixed=TRUE)
  lineRpim <- contents_aimless[gRpim]
  tmpsplit <- strsplit(lineRpim," ")
  idx <- which(nchar(tmpsplit[[1]]) != 0)
  Rpim_values <- as.numeric(tmpsplit[[1]][idx[6:8]]) 
  gCompleteness <- grep("Completeness  ",contents_aimless,fixed=TRUE)
  lineCompleteness <- contents_aimless[gCompleteness]
  tmpsplit <- strsplit(lineCompleteness," ")
  idx <- which(nchar(tmpsplit[[1]]) != 0)
  Completeness_values <- as.numeric(tmpsplit[[1]][idx[2:4]]) 
  stmp <- sprintf("Rmeas, Rpim and Completeness for this cycle: %10.3f, %10.3f, %10.1f\n",Rmeas_values[1],Rpim_values[1],Completeness_values[1])
  cat(stmp)
 
  # Rmerge
  gRmerge <- grep("N   Batch    Mn(I)   RMSdev  I/rms  Rmerge",contents_aimless,fixed=TRUE)
  istart <- gRmerge[1]+1
  iend <- gRmerge[2]-2
  line_Rmerge <- contents_aimless[istart:iend]
  #tbl <- get_vectors_from_table(line_Rmerge,c(2,6,9,13))
  tbl <- get_vectors_from_table(line_Rmerge,c(2,6,9))

  # This bit is just to write out values to be loaded in an interactive R session 
  # for testing 
  #fff <- "./ddd.dat"
  #stmp <- sprintf("%6d %10.3f %8.3f %8.3f\n",tbl[1,1],tbl[1,2],tbl[1,3],tbl[1,4])
  #cat(stmp,file=fff)
  #for (i in 2:length(tbl[,1]))
  #{
  # stmp <- sprintf("%6d %10.3f %8.3f %8.3f\n",tbl[i,1],tbl[i,2],tbl[i,3],tbl[i,4])
  # cat(stmp,file=fff,append=TRUE)
  #}

  # Store average Rmerge and Images per Run
  list_Rmerges <- list()
  list_Images <- list()

  # First exclusion plan
  lista <- pruning_plan(tbl,completeness,cutfraction,0)
  ini_ima <- lista$ini_ima
  fin_ima <- lista$fin_ima
  n_im_elim <- lista$n_im_elim
  neq <- lista$neq
  isel <- lista$isel
  list_Rmerges <- c(list_Rmerges,list(lista$mean_Rmerge))
  list_Images <- c(list_Images,list(lista$Images))
  rm(lista)
  if (neq[isel] == 0)
  {
   stmp <- "=> No images can be further eliminated without either decreasing target completeness or cancelling a whole dataset.\n"
   stmp <- paste(stmp,"  Cycles will be skipped.\n")
   cat(stmp)
  }
  if (neq[isel] > 0)
  {
   cat(paste("=> ",n_im_elim," images can be further eliminated without affecting completeness.\n"))
   cat(paste("   ",neq[isel]," images will be eliminated from run ",isel,".\n"))
   cat("\n")
 
   # Cycles
   icyc <- 0
   repeat
   {
    icyc <- icyc+1
    if (icyc > maxcycle) 
    {
     stmp <- "=> Max. number of cycles reached.\n"
     stmp <- paste(stmp,"Cycles will be stopped\n")
     cat(stmp)
     break
    }
    cat(paste("### Cycle",icyc,"###\n"))
 
    # Load fixed aimless keywords
    gkeys1 <- grep("Input command lines",contents_aimless,fixed=TRUE)
    gkeys2 <- grep("End of input",contents_aimless,fixed=TRUE)
    aimless_keys <- contents_aimless[(gkeys1+1):(gkeys2-1)]
    idx <- which(aimless_keys == "END")
    if (length(idx) != 0) aimless_keys <- aimless_keys[-idx]

    # Additional line for aimless keywords to exclude images
    stmp <- paste("EXCLUDE BATCH",ini_ima,"TO",fin_ima)
    aim_keys <- aimless_keys
    aim_keys <- c(aim_keys,stmp,"END")
    tmp <- simple_merge_datasets(suffix,aim_keys,rwin=rwin)                 ### AIMLESS EXECUTION!!!
    stored_results <- c(stored_results,list(tmp[[1]]))
    stored_aimkeys <- c(stored_aimkeys,list(aim_keys))
    if (is.na(tmp[[1]]$Rmeas[1])) 
    {
     stmp <- "AIMLESS cannot converge to finite results.\n"
     cat(stmp)
     stmp <- "Cycles will be stopped.\n"
     cat(stmp)
     break 
    } 

    # Read AIMLESS log for next cycle and store some info
    log_aimless <- paste(suffix[1],paste("aimless_",suffix[2],".log",sep=""),sep="/")
    contents_aimless <- scan(log_aimless,what="character",sep="\n",quiet=TRUE)
 
    # Rpim and Completeness for this cycle
    gRmeas <- grep("Rmeas (all",contents_aimless,fixed=TRUE)
    lineRmeas <- contents_aimless[gRmeas]
    tmpsplit <- strsplit(lineRmeas," ")
    idx <- which(nchar(tmpsplit[[1]]) != 0)
    Rmeas_values <- as.numeric(tmpsplit[[1]][idx[6:8]]) 
    gRpim <- grep("Rpim (all",contents_aimless,fixed=TRUE)
    lineRpim <- contents_aimless[gRpim]
    tmpsplit <- strsplit(lineRpim," ")
    idx <- which(nchar(tmpsplit[[1]]) != 0)
    Rpim_values <- as.numeric(tmpsplit[[1]][idx[6:8]]) 
    gCompleteness <- grep("Completeness  ",contents_aimless,fixed=TRUE)
    lineCompleteness <- contents_aimless[gCompleteness]
    tmpsplit <- strsplit(lineCompleteness," ")
    idx <- which(nchar(tmpsplit[[1]]) != 0)
    Completeness_values <- as.numeric(tmpsplit[[1]][idx[2:4]]) 
    stmp <- sprintf("Rmeas, Rpim and Completeness for this cycle: %10.3f, %10.3f, %10.1f\n",Rmeas_values[1],Rpim_values[1],Completeness_values[1])
    cat(stmp)

    # If completeness is below the input threshold stop and delete results for this cycle
    if (Completeness_values[1] < completeness)
    {
     stored_results[[length(stored_results)]] <- NULL
     stored_aimkeys[[length(stored_aimkeys)]] <- NULL
     stmp <- "=> No images can be further eliminated without either decreasing target completeness or cancelling a whole dataset.\n"
     stmp <- paste(stmp,"Cycles will be stopped\n")
     cat(stmp)
     break
    }

    # Rmerge and cumulative completeness
    gRmerge <- grep("N   Batch    Mn(I)   RMSdev  I/rms  Rmerge",contents_aimless,fixed=TRUE)
    istart <- gRmerge[1]+1
    iend <- gRmerge[2]-2
    line_Rmerge <- contents_aimless[istart:iend]
    #tbl <- get_vectors_from_table(line_Rmerge,c(2,6,9,13))
    tbl <- get_vectors_from_table(line_Rmerge,c(2,6,9))
  
    # This bit is just to write out values to be loaded in an interactive R session 
    # for testing 
    #fff <- "./ddd.dat"
    #stmp <- sprintf("%6d %10.3f %8.3f %8.3f\n",tbl[1,1],tbl[1,2],tbl[1,3],tbl[1,4])
    #cat(stmp,file=fff)
    #for (i in 2:length(tbl[,1]))
    #{
    # stmp <- sprintf("%6d %10.3f %8.3f %8.3f\n",tbl[i,1],tbl[i,2],tbl[i,3],tbl[i,4])
    # cat(stmp,file=fff,append=TRUE)
    #}
 
    # Exclusion plan for next cycle
    lista <- pruning_plan(tbl,completeness,cutfraction,icyc)
    ini_ima <- lista$ini_ima
    fin_ima <- lista$fin_ima
    n_im_elim <- lista$n_im_elim
    neq <- lista$neq
    isel <- lista$isel
    list_Rmerges <- c(list_Rmerges,list(lista$mean_Rmerge))
    list_Images <- c(list_Images,list(lista$Images))
    rm(lista)
    if (neq[isel] == 0)
    {
     stmp <- "=> No images can be further eliminated without either decreasing target completeness or cancelling a whole dataset.\n"
     stmp <- paste(stmp,"Cycles will be stopped\n")
     cat(stmp)
     break
    }
    if (neq[isel] > 0)
    {
     cat(paste("=> ",n_im_elim," images can be further eliminated at the next cycle without affecting completeness.\n"))
     cat(paste("   ",neq[isel]," images will be eliminated from run ",isel,".\n"))
     cat("\n")
    }
   }
  }

  # Decide which cycle has been the best (for now in terms of Rpim) and re-run AIMLESS
  # job to store related results
  Rpims <- c()
  tCmpl <- c()
  for (icyc in 1:length(stored_results))
  {
   Rpims <- c(Rpims,stored_results[[icyc]]$Rpim[1])
   tCmpl <- c(tCmpl,stored_results[[icyc]]$Completeness[1])
  }
  Rpims[Rpims < 0] <- NA
  idx <- which(Rpims == min(Rpims,na.rm=TRUE))
  if (length(idx > 1))                          # All this (and tCmpl) is because, accidentally, you can have two minimum Rpims
  {
   jdx <- which(tCmpl[idx] == max(tCmpl[idx]))
   idx <- idx[jdx[1]]
  }
  cat("\n")
  cat(paste("Best results have been produced at cycle",(idx-1),". Re-running that cycle ...\n"))
  if (idx == 1)
  {
   aimless_keys <- stored_aimkeys[[idx]]
   tmp <- merge_datasets("FINAL_list_of_files.dat",selection=cln,suffix,pointless_keys,aimless_keys,
                         resomax=resosel,nref=idxref_char,rwin=rwin)
  }
  if (idx > 1)
  {
   aim_keys <- stored_aimkeys[[idx]]
   tmp <- simple_merge_datasets(suffix,aim_keys,rwin=rwin)
  }
 }
 if (is.na(tmp[[1]]$Rmeas[1]))
 {
  cat("\n")
  stmp <- "Starting combination cannot be scaled by AIMLESS (see AIMLESS log). Review keywords and run BLEND again.\n"
  cat(stmp)
 }

 # Prepare and output table with Images and average Rmerge per run at every cycle
 tmp_n <- length(stored_aimkeys)
 tmp_m <- length(list_Images)
 ncycs_final <- min(tmp_n,tmp_m)
 nruns_final <- length(list_Images[[1]])
 mat_Rmerges_final <- matrix(ncol=nruns_final,nrow=ncycs_final)
 mat_Images_final <- matrix(ncol=nruns_final,nrow=ncycs_final)
 for (i in 1:ncycs_final)
 {
  for (j in 1:nruns_final)
  {
   mat_Rmerges_final[i,j] <- list_Rmerges[[i]][j] 
   mat_Images_final[i,j] <- list_Images[[i]][j] 
  }
 }
 cat("\n")
 stmp <- "$TABLE: Image number per cycle for all datasets :\n"
 cat(stmp)
 stmp <- "$GRAPHS:    Number of images per cycle         : N : 1"
 for (j in 2:(nruns_final+1))
 {
  stmp <- paste(stmp,sprintf(", %3d",j),sep="")
 }
 stmp <- paste(stmp," :\n",sep="")
 cat(stmp)
 cat("$$\n")
 stmp <- "  Cycle"
 for (j in 1:nruns_final)
 {
  stmp <- paste(stmp,sprintf("    Dataset_%03d",cln[j]),sep="")
 }
 stmp <- paste(stmp," $$\n",sep="")
 cat(stmp)
 cat("$$\n")
 for (i in 1:ncycs_final)
 {
  stmp <- sprintf("%7d",i)
  for (j in 1:nruns_final)
  {
   stmp <- paste(stmp,sprintf(" %15d",mat_Images_final[i,j]),sep="")
  }
  stmp <- paste(stmp,"\n",sep="")
  cat(stmp)
 }
 cat("$$\n")
 cat("\n")
 stmp <- "$TABLE: Average Rmerge per cycle for all datasets :\n"
 cat(stmp)
 stmp <- "$GRAPHS:    Average Rmerge per cycle         : N : 1"
 for (j in 2:(nruns_final+1))
 {
  stmp <- paste(stmp,sprintf(", %3d",j),sep="")
 }
 stmp <- paste(stmp," :\n",sep="")
 cat(stmp)
 cat("$$\n")
 stmp <- "  Cycle"
 for (j in 1:nruns_final)
 {
  stmp <- paste(stmp,sprintf("    Dataset_%03d",cln[j]),sep="")
 }
 stmp <- paste(stmp," $$\n",sep="")
 cat(stmp)
 cat("$$\n")
 for (i in 1:ncycs_final)
 {
  stmp <- sprintf("%7d",i)
  for (j in 1:nruns_final)
  {
   stmp <- paste(stmp,sprintf("%15.3f",mat_Rmerges_final[i,j]),sep="")
  }
  stmp <- paste(stmp,"\n",sep="")
  cat(stmp)
 }
 cat("$$\n")
 cat("\n")
}
if (combination_type == 2)
{
 # Group of datasets is constantly updated
 current_cln <- cln

 # Create two lists to store aimless keywords and merging statistics.
 stored_results <- list()
 stored_selections <- list()

 # Initial run of AIMLESS
 cat("### Cycle 0 ###\n")
 cat("No filtering for the initial cycle.\n")
 tmp <- merge_datasets("FINAL_list_of_files.dat",selection=current_cln,suffix,pointless_keys,aimless_keys,
                       resomax=resosel,nref=idxref_char,rwin=rwin)

 # Get table of average Rmerges per run
 if (!is.na(tmp[[1]]$Rmeas[1]))
 {
  stored_results <- c(stored_results,list(tmp[[1]]))
  stored_selections <- c(stored_selections,list(current_cln))

  # Read AIMLESS log and store some info
  log_aimless <- paste(suffix[1],paste("aimless_",suffix[2],".log",sep=""),sep="/")
  contents_aimless <- scan(log_aimless,what="character",sep="\n",quiet=TRUE)
 
  # Rpim and Completeness for this (starting) cycle
  gRmeas <- grep("Rmeas (all",contents_aimless,fixed=TRUE)
  lineRmeas <- contents_aimless[gRmeas]
  tmpsplit <- strsplit(lineRmeas," ")
  idx <- which(nchar(tmpsplit[[1]]) != 0)
  Rmeas_values <- as.numeric(tmpsplit[[1]][idx[6:8]]) 
  gRpim <- grep("Rpim (all",contents_aimless,fixed=TRUE)
  lineRpim <- contents_aimless[gRpim]
  tmpsplit <- strsplit(lineRpim," ")
  idx <- which(nchar(tmpsplit[[1]]) != 0)
  Rpim_values <- as.numeric(tmpsplit[[1]][idx[6:8]])
  gCompleteness <- grep("Completeness  ",contents_aimless,fixed=TRUE)
  lineCompleteness <- contents_aimless[gCompleteness]
  tmpsplit <- strsplit(lineCompleteness," ")
  idx <- which(nchar(tmpsplit[[1]]) != 0)
  Completeness_values <- as.numeric(tmpsplit[[1]][idx[2:4]])
  stmp <- sprintf("Rmeas, Rpim and Completeness for this cycle: %10.3f, %10.3f, %10.1f\n",Rmeas_values[1],Rpim_values[1],Completeness_values[1])
  cat(stmp)

  # Rmerge
  gRmerge <- grep("N   Batch    Mn(I)   RMSdev  I/rms  Rmerge",contents_aimless,fixed=TRUE)
  istart <- gRmerge[1]+1
  iend <- gRmerge[2]-2
  line_Rmerge <- contents_aimless[istart:iend]
  #tbl <- get_vectors_from_table(line_Rmerge,c(2,6,9,13))
  tbl <- get_vectors_from_table(line_Rmerge,c(2,6,9))

  # Store average Rmerge and Images per Run
  list_Rmerges <- list()
  list_Images <- list()

  # First exclusion plan
  lista <- pruning_plan(tbl,100,1,0)
  mean_Rmerge <- abs(lista$mean_Rmerge)
  list_Rmerges <- c(list_Rmerges,list(mean_Rmerge))
  list_Images <- c(list_Images,list(lista$Images))
  rm(lista)
 
  # Cycle through decreasing datasets
  icyc <- 0
  repeat
  {
   icyc <- icyc+1
   if (icyc > maxcycle)
   {
    stmp <- "=> Max. number of cycles reached.\n"
    stmp <- paste(stmp,"Cycles will be stopped\n")
    cat(stmp)
    break
   }
   if (length(mean_Rmerge) == 1)
   {
    stmp <- "=> No more datasets left to be processed. Filtering completed\n"
    cat(stmp)
    break
   }
   cat("\n")
   cat(paste("### Cycle",icyc,"###\n"))

   # Determine which dataset produced by the previous AIMLESS run has worst Rmerge
   idxrun <- which(mean_Rmerge == max(mean_Rmerge,na.rm=TRUE))[1]
   stmp <- sprintf("Filtering out dataset %d.\n",current_cln[idxrun])
   cat(stmp)

   # Eliminate dataset from combination and re-run POINTLESS+AIMLESS
   current_cln <- current_cln[-idxrun]
   tmp <- merge_datasets("FINAL_list_of_files.dat",selection=current_cln,suffix,pointless_keys,aimless_keys,
                         resomax=resosel,nref=idxref_char,rwin=rwin)
   stored_results <- c(stored_results,list(tmp[[1]]))
   stored_selections <- c(stored_selections,list(current_cln))
   if (is.na(tmp[[1]]$Rmeas[1]))
   {
    stmp <- "AIMLESS cannot converge to finite results.\n"
    cat(stmp)
    stmp <- "Cycles will be stopped.\n"
    cat(stmp)
    break
   }
 
   # Read AIMLESS log for next cycle and store some info
   log_aimless <- paste(suffix[1],paste("aimless_",suffix[2],".log",sep=""),sep="/")
   contents_aimless <- scan(log_aimless,what="character",sep="\n",quiet=TRUE)

   # Rpim and Completeness for this cycle
   gRmeas <- grep("Rmeas (all",contents_aimless,fixed=TRUE)
   lineRmeas <- contents_aimless[gRmeas]
   tmpsplit <- strsplit(lineRmeas," ")
   idx <- which(nchar(tmpsplit[[1]]) != 0)
   Rmeas_values <- as.numeric(tmpsplit[[1]][idx[6:8]]) 
   gRpim <- grep("Rpim (all",contents_aimless,fixed=TRUE)
   lineRpim <- contents_aimless[gRpim]
   tmpsplit <- strsplit(lineRpim," ")
   idx <- which(nchar(tmpsplit[[1]]) != 0)
   Rpim_values <- as.numeric(tmpsplit[[1]][idx[6:8]])
   gCompleteness <- grep("Completeness  ",contents_aimless,fixed=TRUE)
   lineCompleteness <- contents_aimless[gCompleteness]
   tmpsplit <- strsplit(lineCompleteness," ")
   idx <- which(nchar(tmpsplit[[1]]) != 0)
   Completeness_values <- as.numeric(tmpsplit[[1]][idx[2:4]])
   stmp <- sprintf("Rmeas, Rpim and Completeness for this cycle: %10.3f, %10.3f, %10.1f\n",Rmeas_values[1],Rpim_values[1],Completeness_values[1])
   cat(stmp)

   # If completeness is below the input threshold stop and delete results for this cycle
   if (Completeness_values[1] < completeness)
   {
    stored_results[[length(stored_results)]] <- NULL
    stored_selections[[length(stored_selections)]] <- NULL
    stmp <- "=> This datasets cannot be eliminated without decreasing target completeness.\n"
    stmp <- paste(stmp,"Cycles will be stopped\n")
    cat(stmp)
    break
   }

   # Rmerge
   gRmerge <- grep("N   Batch    Mn(I)   RMSdev  I/rms  Rmerge",contents_aimless,fixed=TRUE)
   istart <- gRmerge[1]+1
   iend <- gRmerge[2]-2
   line_Rmerge <- contents_aimless[istart:iend]
   #tbl <- get_vectors_from_table(line_Rmerge,c(2,6,9,13))
   tbl <- get_vectors_from_table(line_Rmerge,c(2,6,9))

   # Exclusion plan for next cycle
   lista <- pruning_plan(tbl,100,1,icyc)
   mean_Rmerge <- abs(lista$mean_Rmerge)
   list_Rmerges <- c(list_Rmerges,list(mean_Rmerge))
   list_Images <- c(list_Images,list(lista$Images))
   rm(lista)
  }

  # Collect all results and decide which one to select
  Rpims <- c()
  tCmpl <- c()
  for (icyc in 1:length(stored_results))
  {
   Rpims <- c(Rpims,stored_results[[icyc]]$Rpim[1])
   tCmpl <- c(tCmpl,stored_results[[icyc]]$Completeness[1])
  }
  Rpims[Rpims < 0] <- NA
  idx <- which(Rpims == min(Rpims,na.rm=TRUE))
  if (length(idx > 1))                          # All this (and tCmpl) is because, accidentally, you can have two minimum Rpims
  {
   jdx <- which(tCmpl[idx] == max(tCmpl[idx]))
   idx <- idx[jdx[1]]
  }
  cat("\n")
  cat(paste("Best results have been produced at cycle",(idx-1),". Re-running that cycle ...\n"))
  current_cln <- stored_selections[[idx]]
  tmp <- merge_datasets("FINAL_list_of_files.dat",selection=current_cln,suffix,pointless_keys,aimless_keys,
                        resomax=resosel,nref=idxref_char,rwin=rwin)

  # Output table with Rpim and Rmeas for each cycle
  ncycs <- length(stored_results)
  cat("\n")
  cat("$TABLE: Rmeas and Rpim for all cycles :\n")
  cat("$GRAPHS:    Overall Rmeas and Rpim        : N : 1, 2, 3 :\n")
  cat("$$\n")
  cat("     Cycle           Rmeas            Rpim  $$\n")
  cat("$$\n")
  for (icyc in 1:ncycs)
  {
   stmp <- sprintf("%10d %15.3f %15.3f\n",(icyc-1),stored_results[[icyc]]$Rmeas[1],stored_results[[icyc]]$Rpim[1])
   cat(stmp)
  }
  cat("$$\n")
  cat("\n")

  # Before moving on to final part of the program, copy chosen cln to cln
  cln <- current_cln 
 }

 # Terminate if no results come from initial scaling
 if (is.na(tmp[[1]]$Rmeas[1]))
 {
  cat("\n")
  stmp <- "Starting combination cannot be scaled by AIMLESS (see AIMLESS log). Review keywords and run BLEND again.\n"
  cat(stmp)
 }
}

# Carry on with displaying results
cat("\n")
cat(" Statistics for this group:\n")

# Change row names for display purpose
rownames(tmp[[1]]) <- c("Overall","InnerShell","OuterShell")
print(tmp[[1]])
if (is.na(tmp[[1]][1,1])) cat(paste("WARNING! No result could be produced for this group, due to a problem with either POINTLESS or AIMLESS.\n",
                                       "         You can try and explore why in directory <combined_files> (see documentation).\n\n",sep=""))

# Add to what's already included in "combined_files" directory
cat("\n")
if (file.exists(outdir))
{
 # If object groups_info exists, create file GROUPS.info and dump existing stuff
 if (exists("groups_info"))
 {
  mtmp <- paste(outdir,"/GROUPS.info",sep="")
  if (file.exists(mtmp)) emptyc <- file.remove(mtmp)
  linea <- paste(groups_info[1],"\n",sep="")
  cat(linea,file=mtmp)
  for (i in 2:length(groups_info))
  {
   linea <- paste(groups_info[i],"\n",sep="")
   cat(linea,file=mtmp,append=TRUE)
  }
 }

 # If object merging_statistics_info exists, create file MERGING_STATISTICS.info and dump existing stuff
 if (exists("merging_statistics_info"))
 {
  mtmp <- paste(outdir,"/MERGING_STATISTICS.info",sep="")
  if (file.exists(mtmp)) emptyc <- file.remove(mtmp)
  linea <- paste(merging_statistics_info[1],"\n",sep="")
  cat(linea,file=mtmp)
  for (i in 2:length(merging_statistics_info))
  {
   linea <- paste(merging_statistics_info[i],"\n",sep="")
   cat(linea,file=mtmp,append=TRUE)
  }
 }

 # Now add this result to existing files, or write new files if first result
 if (!exists("groups_info"))
 {
  mtmp <- paste(outdir,"/GROUPS.info",sep="")
  linea <- sprintf("********* Group %d *********\n",(length(gidx)+1))
  cat(linea,file=mtmp,append=TRUE)
  for (k in cln)
  {
   jcn <- which(maindf$cn == k)
   linea <- sprintf("%-100s %5.0f   %5.0f\n",filenames[as.integer(k),],maindf[jcn,"cn"],maindf[jcn,"batch"])
   cat(linea,file=mtmp,append=T)
  }
 }
 if (exists("groups_info"))
 {
  mtmp <- paste(outdir,"/GROUPS.info",sep="")
  linea <- sprintf("********* Group %d *********\n",(length(gidx)+1))
  cat(linea,file=mtmp,append=T)
  for (k in cln)
  {
   jcn <- which(maindf$cn == k)
   linea <- sprintf("%-100s %5.0f   %5.0f\n",filenames[as.integer(k),],maindf[jcn,"cn"],maindf[jcn,"batch"])
   cat(linea,file=mtmp,append=T)
  }
 }

 # Extract CC1/2 and Mn2 and add them to tmp (other statistics)
 if (!is.null(tmp[[2]]) & length(tmp[[2]]) > 0)
 {
  gCC12 <- grep("from half-dataset correlation CC(1/2)",tmp[[2]],fixed=TRUE)
  lineCC12 <- tmp[[2]][gCC12][1]
  stmp <- strsplit(lineCC12,">")[[1]][2]
  stmp <- strsplit(stmp,"=")
  stmp <- gsub("\\s","",stmp[[1]][2])
  lstmp <- strsplit(stmp,"A")
  CC12 <- as.numeric(lstmp[[1]][1])
 }
 if (is.null(tmp[[2]]) | length(tmp[[2]]) == 0) CC12 <- NA
 
 # Extract Mn2
 if (!is.null(tmp[[2]]) & length(tmp[[2]]) > 0)
 {
  gMn2 <- grep("from Mn(I/sd) >  1.50:",tmp[[2]],fixed=TRUE)
  lineMn2 <- tmp[[2]][gMn2][1]
  stmp <- strsplit(lineMn2,"=")[[1]][2]
  stmp <- gsub("\\s","",stmp)
  Mn2 <- as.numeric(substr(stmp,1,(nchar(stmp)-1)))
 }
 if (is.null(tmp[[2]]) | length(tmp[[2]]) == 0) Mn2 <- NA
 
 extra_tmp <- data.frame(CC12=CC12,Mn2=Mn2)

 logdframe <- data.frame()
 if (!exists("merging_statistics_info"))
 {
  mtmp <- paste(outdir,"/MERGING_STATISTICS.info",sep="")
  linea <- sprintf("######################################################################################################################\n")
  cat(linea,file=mtmp) 
  linea <- sprintf("######################################################################################################################\n")
  cat(linea,file=mtmp,append=TRUE) 
  linea <- sprintf("###################################### MERGING STATISTICS FOR ALL INPUT GROUPS #######################################\n")
  cat(linea,file=mtmp,append=TRUE) 
  linea <- sprintf("######################################################################################################################\n")
  cat(linea,file=mtmp,append=TRUE) 
  linea <- sprintf("######################################################################################################################\n")
  cat(linea,file=mtmp,append=TRUE) 
  linea <- "                            \n"
  cat(linea,file=mtmp,append=TRUE)
  linea <- sprintf("   Group                     Comple-    Multi-     Reso     Reso      Reso\n")
  cat(linea,file=mtmp,append=T)
  linea <- sprintf("  Number    Rmeas     Rpim   teness    plicity    CC1/2   Mn(I/sd)     Max\n")
  cat(linea,file=mtmp,append=T)
  linea <- sprintf("     %3d  %7.3f  %7.3f   %6.2f     %6.2f  %7.2f   %7.2f  %7.2f\n",
                   (length(gidx)+1),tmp[[1]][1,1],tmp[[1]][1,2],tmp[[1]][1,3],tmp[[1]][1,4],extra_tmp[1,1],extra_tmp[1,2],tmp[[1]][1,6])
  cat(linea,file=mtmp,append=TRUE)
  logdframe <- rbind(logdframe,cbind(1,tmp[[1]][1,c(1,2,3,4)],extra_tmp[1,],tmp[[1]][1,6]))
 }
 if (exists("merging_statistics_info"))
 {
  mtmp <- paste(outdir,"/MERGING_STATISTICS.info",sep="")
  linea <- sprintf("     %3d  %7.3f  %7.3f   %6.2f     %6.2f  %7.2f   %7.2f  %7.2f\n",
                  (length(gidx)+1),tmp[[1]][1,1],tmp[[1]][1,2],tmp[[1]][1,3],tmp[[1]][1,4],extra_tmp[1,1],extra_tmp[1,2],tmp[[1]][1,6])
  cat(linea,file=mtmp,append=TRUE)
  logdframe <- rbind(logdframe,read.table(mtmp,skip=8))
 }
}

# Output for logview
tmpdframe <- logdframe
if (length(logdframe[,1]) > 1) logdframe <- tmpdframe[order(-tmpdframe[,4],tmpdframe[,2],na.last=TRUE),]
if (length(logdframe[,1]) == 1) rownames(logdframe) <- "1"                                                                   # Fix logdframe temporarily
linea <- sprintf("$TABLE: Overall merging statistics and completeness :\n")                                                  # For logview
cat(linea)                                                                                                                   # For logview
linea <- sprintf("$GRAPHS:        Overall merging statistics   : N : 1, 3, 4 :\n")                                           # For logview
cat(linea)                                                                                                                   # For logview
linea <- sprintf("       :        Overall completeness         : N : 1, 5 :\n")                                              # For logview
cat(linea)                                                                                                                   # For logview
linea <- sprintf("       :        Resolutions (CC1/2, Mn2, Max): N : 1, 7, 8, 9 :\n")                                              # For logview
cat(linea)                                                                                                                   # For logview
linea <- sprintf("$$\n")                                                                                                     # For logview
cat(linea)                                                                                                                   # For logview
linea <- sprintf(" Index    Group    Rmeas     Rpim   Compl.   Multip.   Res_CC1/2    Res_Mn(I/sd)   Res_Max  $$\n")                    # For logview
cat(linea)
linea <- sprintf("$$\n")                                                                                                     # For logview
cat(linea)                                                                                                                   # For logview
for (i in 1:length(logdframe[,1]))
{
 j <- as.integer(rownames(logdframe)[i])
 linea2 <- sprintf("   %3d      %3d  %7.3f  %7.3f  %6.2f    %6.2f      %7.2f        %7.2f    %7.2f\n",
                   i,logdframe[i,1],logdframe[i,2],logdframe[i,3],logdframe[i,4],logdframe[i,5],logdframe[i,6],logdframe[i,7],logdframe[i,8])
 cat(linea2)                                                                                                                 # For logview
}
linea <- sprintf("$$\n")                                                                                                     # For logview
cat(linea)                                                                                                                   # For logview
linea <- sprintf("\n")                                                                                                       # For logview
cat(linea)

# Exit without saving
q(save = "no", status = 0, runLast = FALSE)
