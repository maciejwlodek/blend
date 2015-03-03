#***********************************************************************************************************
#***********************************************************************************************************
#********* blend2.R
#*********
#********* Copyright (C) 2014 Diamond Light Source & Imperial College London
#*********
#********* Authors: James Foadi & Gwyndaf Evans
#*********
#********* This code is distributed under the BSD license, a copy of which is
#********* included in the root directory of this package.
#************************************************************************************************************
#************************************************************************************************************
# Second R module of program BLEND (later to be replaced by C++ code)
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
   linea <- sprintf("EXCLUDE FILE %d BATCH %d TO %d\n",imtz,inibatch,finbatch)
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
  for (linea in exemerge)
  {
   linea <- paste(linea,"\n",sep="")
   cat(linea,file=log_file,append=TRUE)
  }

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
   idx <- grep("reference",dircontents,fixed=TRUE)
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
  for (linea in exeaimless)
  {
   linea <- paste(linea,"\n",sep="")
   cat(linea,file=log_file,append=TRUE)
  }

  # Remove files produced to merge mtz's
  dircontents <- list.files("./")
  idx <- grep("mtzdump",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (i in idx) file.remove(dircontents[i])
  idx <- grep("pointless",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (i in idx) file.remove(dircontents[i])
  idx <- grep("reference",dircontents,fixed=TRUE)
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
 else
 {
  Mstats <- data.frame(Rmeas=c(NA,NA,NA),Rpim=c(NA,NA,NA),Completeness=c(NA,NA,NA),Multiplicity=c(NA,NA,NA),
                       LowRes=c(NA,NA,NA),HighRes=c(NA,NA,NA))

  # Rename POINTLESS log
  log_file <- paste(suffix[1],paste("pointless_",suffix[2],".log",sep=""),sep="/")
  for (linea in exemerge)
  {
   linea <- paste(linea,"\n",sep="")
   cat(linea,file=log_file,append=TRUE)
  }
  
  return(list(Mstats,exemerge))
 }
}



######################################################################################################################
######################################################################################################################

# Main

# Find out if running on Windows
ostuff <- Sys.info()
rwin = FALSE
if (ostuff["sysname"] == "Windows") rwin = TRUE

# Load content of previous R run
load("BLEND.RData",.GlobalEnv)

# To avoid warning messages set warn to a negative value 
options(warn = -1)

# Retrieve value from command line
args <- commandArgs(trailingOnly=TRUE)
tmp <- as.numeric(args)
level_type <- as.integer(tmp[1])
dlevel_top <- tmp[2]
dlevel_bottom <- tmp[3]

# List of wanted nodes (from input)
if (level_type == 0)
{
 if (dlevel_bottom < 0) idx <- which(npar.hc_ward$height < dlevel_top)
 if (dlevel_bottom > 0) idx <- which(npar.hc_ward$height > dlevel_bottom & npar.hc_ward$height < dlevel_top)
}
if (level_type == 1)
{
 if (dlevel_bottom < 0) idx <- which(LCV_values < dlevel_top)
 if (dlevel_bottom > 0) idx <- which(LCV_values > dlevel_bottom & LCV_values < dlevel_top)
}
if (level_type == 2)
{
 if (dlevel_bottom < 0) idx <- which(aLCV_values < dlevel_top)
 if (dlevel_bottom > 0) idx <- which(aLCV_values > dlevel_bottom & aLCV_values < dlevel_top)
}

# If no merging nodes are included between the selected levels, stop with a warning
if (length(idx) == 0)
{
 if (dlevel_bottom < 0) dlevel_bottom <- 0.0
 cat("\n")
 cat(paste("WARNING! Your selection range (",dlevel_bottom," - ",dlevel_top,") includes no merging nodes. Try a different range.\n",sep=""))
 cat("\n")
}
if (length(idx) > 0)
{
 # Extract keywords files for POINTLESS and AIMLESS
 pointless_keys <- c()
 aimless_keys <- c()
 if (file.exists("BLEND_KEYWORDS.dat"))
 {
  contents <- scan("BLEND_KEYWORDS.dat", what="character",sep="\n",quiet=TRUE)

  # POINTLESS
  idxs <- grep("POINTLESS KEYWORDS",contents,fixed=TRUE)
  idxe <- grep("AIMLESS KEYWORDS",contents,fixed=TRUE)
  if ((idxe-idxs) > 1) for (line in contents[(idxs+1):(idxe-1)]) pointless_keys <- c(pointless_keys,line)

  # Reference dataset (from BLEND KEYWORDS section)
  tmp <- grep("DATAREF",contents,fixed=TRUE)
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

  # AIMLESS
  idxs <- idxe
  idxe <- length(contents)
  if ((idxe-idxs) > 0) for (line in contents[(idxs+1):idxe]) aimless_keys <- c(aimless_keys,line)
 }
 
 # Create directory to store clusters file. If directory exists, delete previous content.
 outdir <- "./merged_files"
 if (!file.exists(outdir)) dir.create(outdir)
 if (file.exists(outdir))
 {
  #dircontents <- system(paste("ls ",outdir,sep=""),intern=TRUE)
  dircontents <- list.files(outdir)
  for (a in dircontents)
  {
   line <- paste(outdir,a,sep="/")
   emtyc <- file.remove(line)
  }
 }
 
 # Ascii file with list of merging datasets
 clusters_list_file=paste(outdir,"CLUSTERS.info",sep="/")
 if (file.exists(clusters_list_file)) emptyc <- file.remove(clusters_list_file)

 # Loop over all selected merging nodes and produce scaled datasets
 mergingStatistics <- data.frame()
 for (i in 1:length(idx))
 {
  j <- idx[i]
 
  # Produce text related to specific cluster in CLUSTERs.info file
  #cln <- groups[[1]][[j]]
  tmpcln <- groups[[1]][[j]]
  cln <- sort(tmpcln)
  messaggio <- paste("********* Cluster ",j,", composed of datasets ",sep="")
  for (k in cln) messaggio <- paste(messaggio,k," ",sep="")
  messaggio <- paste(messaggio,"*********\n",sep="")
  cat(messaggio)
  if (length(cln) > 1)
  {
   linea <- sprintf("********* Cluster %d *********\n",j)
   cat(linea,file=clusters_list_file,append=T)
   for (k in cln)
   {
    jcn <- which(maindf$cn == k)
    linea <- sprintf("%-100s %5.0f   %5.0f\n",filenames[as.integer(k),],maindf[jcn,"cn"],maindf[jcn,"batch"])
    cat(linea,file=clusters_list_file,append=T)
   }
  
   # Scale and merge datasets in this specific cluster
   suffix <- c(outdir,sprintf("%03d",j))
   tmp <- merge_datasets("FINAL_list_of_files.dat",selection=cln,suffix,pointless_keys,aimless_keys,
                         resomin=groups[[2]][[j]][1],resomax=groups[[2]][[j]][2],nref=idxref_char,rwin=rwin)
   cat(" Statistics for this group:\n")

   # Change row names for display purpose
   rownames(tmp[[1]]) <- c("Overall","InnerShell","OuterShell")
   print(tmp[[1]])
   if (is.na(tmp[[1]][1,1])) cat(paste("WARNING! No result could be produced for cluster ",j," due to a problem with either POINTLESS or AIMLESS.\n",
                                       "         You can try and explore why in directory <merged_files> (see documentation).\n\n",sep=""))
 
   # Extract CC1/2
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

   # Extract Reso Mn(I/sd)
   if (!is.null(tmp[[2]]) & length(tmp[[2]]) > 0)
   {
    gMn2 <- grep("from Mn(I/sd) >  1.50:",tmp[[2]],fixed=TRUE)
    lineMn2 <- tmp[[2]][gMn2][1]
    stmp <- strsplit(lineMn2,"=")[[1]][2]
    stmp <- gsub("\\s","",stmp)
    Mn2 <- as.numeric(substr(stmp,1,(nchar(stmp)-1)))
   }
   if (is.null(tmp[[2]]) | length(tmp[[2]]) == 0) Mn2 <- NA
   
   # Collect merging statistics in data frame
   mergingStatistics <- rbind(mergingStatistics,data.frame(tmp[[1]][1,],CC12=CC12,Mn2=Mn2))
  }
 }
 rownames(mergingStatistics) <- idx

 save(mergingStatistics,file="BLEND.RMergingStatistics")

 # Plot Rmeas vs Completeness (in both PNG and POSTSCRIPT formats)
 tmp <- is.na(mergingStatistics$Rmeas)
 ntmp <- sum(tmp)
 if (ntmp < length(mergingStatistics[,1]))
 {
  png("./merged_files/Rmeas_vs_Cmpl.png")
  plot(mergingStatistics$Completeness,mergingStatistics$Rmeas,cex=3,xlab="Completeness",ylab="Rmeas")
  text(mergingStatistics$Completeness,mergingStatistics$Rmeas,labels=rownames(mergingStatistics),cex=1)
  emptyc <- dev.off()
  postscript("./merged_files/Rmeas_vs_Cmpl.ps", width = 10, height = 10, paper = "a4")
  plot(mergingStatistics$Completeness,mergingStatistics$Rmeas,cex=3,xlab="Completeness",ylab="Rmeas")
  text(mergingStatistics$Completeness,mergingStatistics$Rmeas,labels=rownames(mergingStatistics),cex=1)
  emptyc <- dev.off()
 }
 
 # Sort mergingStatistics data frame according to highest completeness and lowest Rmeas
 tmpdframe <- mergingStatistics
 #mergingStatistics <- tmpdframe[order(tmpdframe$Rmeas,tmpdframe$Completeness,na.last=TRUE),]
 mergingStatistics <- tmpdframe[order(-tmpdframe$Completeness,tmpdframe$Rmeas,na.last=TRUE),]

 # Output final table (output to log file as well)
 merging_statistics_file=paste(outdir,"MERGING_STATISTICS.info",sep="/")
 if (file.exists(merging_statistics_file)) emptyc <- file.remove(merging_statistics_file)
 linea <- sprintf("######################################################################################################################\n")
 cat(linea,file=merging_statistics_file,append=T)
 linea <- sprintf("######################################################################################################################\n")
 cat(linea,file=merging_statistics_file,append=T)
 linea <- sprintf("################################## MERGING STATISTICS FOR ALL SELECTED CLUSTERS ######################################\n")
 cat(linea,file=merging_statistics_file,append=T)
 linea <- sprintf("######################################################################################################################\n")
 cat(linea,file=merging_statistics_file,append=T)
 linea <- sprintf("######################################################################################################################\n")
 cat(linea,file=merging_statistics_file,append=T)
 linea <- "                            \n"
 cat(linea,file=merging_statistics_file,append=T)
 linea <- sprintf(" Cluster                     Comple-    Multi-     Reso     Reso      Reso\n")
 cat(linea,file=merging_statistics_file,append=T)
 linea <- sprintf("  Number    Rmeas     Rpim   teness    plicity    CC1/2   Mn(I/sd)     Max\n")
 cat(linea,file=merging_statistics_file,append=T)
 linea <- sprintf("\n")                                                                                                       # For logview
 cat(linea)                                                                                                                   # For logview
 linea <- sprintf("\n")                                                                                                       # For logview
 cat(linea)                                                                                                                   # For logview
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
 linea <- sprintf(" Index  Cluster    Rmeas     Rpim   Compl.   Multip.   Res_CC1/2    Res_Mn(I/sd)   Res_Max  $$\n")                    # For logview
 cat(linea)                                                                                                                   # For logview
 linea <- sprintf("$$\n")                                                                                                     # For logview
 cat(linea)                                                                                                                   # For logview
 for (i in 1:length(mergingStatistics[,1]))
 {
  linea <- sprintf("     %3d  %7.3f  %7.3f   %6.2f     %6.2f  %7.2f   %7.2f  %7.2f\n",
                   as.integer(rownames(mergingStatistics)[i]),
                   mergingStatistics[i,1],mergingStatistics[i,2],mergingStatistics[i,3],mergingStatistics[i,4],
                   mergingStatistics[i,7],mergingStatistics[i,8],mergingStatistics[i,6])
  cat(linea,file=merging_statistics_file,append=T)
  linea2 <- sprintf("   %3d      %3d  %7.3f  %7.3f  %6.2f    %6.2f      %7.2f        %7.2f    %7.2f\n",
                   i,
                   as.integer(rownames(mergingStatistics)[i]),
                   mergingStatistics[i,1],mergingStatistics[i,2],mergingStatistics[i,3],mergingStatistics[i,4],
                   mergingStatistics[i,7],mergingStatistics[i,8],mergingStatistics[i,6])
  cat(linea2)                                                                                                                 # For logview
 }
 linea <- sprintf("$$\n")                                                                                                     # For logview
 cat(linea)                                                                                                                   # For logview
 linea <- sprintf("\n")                                                                                                       # For logview
 cat(linea)                                                                                                                   # For logview

 # Final message for users
 #cat("\n")
 #cat("##################################################################\n")
 #cat("##################################################################\n")
 #cat("##################################################################\n")
 #cat("  The following files have been created in directory merged_files:\n")
 #cat("CLUSTERS.info                         : ascii file including details (file path, serial number and batch used) for each cluster selected;\n")
 #cat("MERGING_STATISTICS.info               : ascii file including a table with merging statistics for each merged and/or scaled file, sorted\n")
 #cat("                                        according to increasing values of Rmeas;\n")
 #cat("Rmeas_vs_Cmpl.png                     : png image displaying Rmeas versus Completeness for all merged nodes\n")
 #cat("merged_001.mtz, merged_002.mtz, etc   : mtz files prior to scaling. One for each cluster;\n")
 #cat("scaled_001.mtz, scaled_002.mtz, etc   : mtz scaled files. One for each cluster;\n")
 #cat("aimless_001.log, aimless_002.log, etc : log files, one for each AIMLESS job.\n")
 #cat("\n")
}

# Exit without saving
q(save = "no", status = 0, runLast = FALSE) 
