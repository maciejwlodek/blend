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
mtz_batch_list <- function(hklin)
{
 # Run mtzdmp and save output
 exemtzdmp <- system(paste("mtzdmp ",hklin,sep=""),intern=TRUE)

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
 
 return(blist)
}

#
# Function to truncate and renumber batches
trunc_and_renum <- function(hklin,inibatch,finbatch,irenum,hklout)
{
 # Build keywords file
 linea <- sprintf("TITLE Rebatch file for multicrystal merging\n")
 cat(linea,file="rebatch_keywords.dat")
 if (inibatch > 0)
 {
  linea <- sprintf("BATCH %d TO %d REJECT\n",inibatch,finbatch)
  cat(linea,file="rebatch_keywords.dat",append=TRUE)
 }
 linea <- sprintf("BATCH ALL START %d\n",irenum)
 cat(linea,file="rebatch_keywords.dat",append=TRUE)
 linea <- sprintf("END\n")
 cat(linea,file="rebatch_keywords.dat",append=TRUE)

 # Run program
 stringa <- paste("rebatch HKLIN ",hklin," HKLOUT ",hklout," < rebatch_keywords.dat")
 exerebatch <- system(stringa,intern=TRUE)

 # Delete keywords file
 if (file.exists("rebatch_keywords.dat")) emptyc <- file.remove("rebatch_keywords.dat")

 return(exerebatch)
}

#
# Alternate indexing with POINTLESS
altidx <- function(hklin,hklref,hklout)
{
 # Pointless keywords file
 linea <- sprintf("NAME PROJECT xxx CRYSTAL yyy DATASET zzz\n")
 cat(linea,file="pointless_keywords.dat")
 linea <- sprintf("HKLIN %s\n",hklin)
 cat(linea,file="pointless_keywords.dat",append=TRUE)
 linea <- sprintf("HKLREF %s\n",hklref)
 cat(linea,file="pointless_keywords.dat",append=TRUE)
 linea <- sprintf("HKLOUT %s\n",hklout)
 cat(linea,file="pointless_keywords.dat",append=TRUE)

 # Run program
 stringa <- "pointless < pointless_keywords.dat"
 exepointless <- system(stringa,intern=TRUE, ignore.stderr = TRUE)

 # Delete keywords file
 if (file.exists("pointless_keywords.dat")) emptyc <- file.remove("pointless_keywords.dat")

 return(exepointless)
}

#
# Run pointless to find Laue group and merge files (after alternate indexing)
merge_mtzs <- function(mtz_list,selection,mtzout,pointless_keys,hklref)
{
 # List of files to remove at the end
 files_to_rm <- c()

 # Cycle through each file in the list and do the following:
 #    1) extract batch information from each mtz file;
 #    2) truncate end of batches according to selection
 #    3) do alternate indexing;
 #    4) glue all mtz together into a single mtz

 #mtzin <- "new_001.mtz"                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 mtzin <- c()
 for (imtz in 1:length(mtz_list))
 {
  # Serial number appended to mtz file
  sfx <- sprintf("_%03d",imtz)

  hklin <- mtz_list[imtz]

  # Extract batch information from each mtz file
  blist <- mtz_batch_list(hklin)

  # Truncate batches
  itmp <- selection[imtz]+1
  if (itmp > blist[length(blist)]) inibatch <- -1
  if (itmp <= blist[length(blist)]) inibatch <- itmp
  finbatch <- blist[length(blist)]
  hklout <- paste("new",sfx,".mtz",sep="")
  files_to_rm <- c(files_to_rm,hklout)
  tmp <- trunc_and_renum(hklin,inibatch,finbatch,1,hklout) 

  # Do alternate indexing
  hklin <- hklout
  hklout <- paste("final",sfx,".mtz",sep="")
  mtzin <- c(mtzin,hklout)
  files_to_rm <- c(files_to_rm,hklout)
  tmp <- altidx(hklin,hklref,hklout)
  
  # Return NULL if space group of reference file is incompatible with space group of any of current files
  exealtidx <- grep("$$ <!--SUMMARY_END-->",tmp,fixed=TRUE)
  if (length(exealtidx) == 0)
  {
   messaggio <- "Reference file is incompatible with some data sets of this cluster.\n"
   cat(messaggio)
   exepointless <- NULL

   # Clean directory from unnecessary files
   for (a in files_to_rm) if (file.exists(a)) file.remove(a)

   return(exepointless)
  }
 }

 # Run Pointless to glue everything together under correct space group
 linea <- sprintf("NAME PROJECT xxx CRYSTAL yyy DATASET zzz\n")
 cat(linea,file="pointless_keywords.dat")
 for (a in mtzin)
 {
  linea <- sprintf("HKLIN %s\n",a)
  cat(linea,file="pointless_keywords.dat",append=TRUE)
 }

 # Additional keywords, introduced by user through BLEND_KEYWORDS.dat file
 for (line in pointless_keys) cat(paste(line,"\n",sep=""),file="pointless_keywords.dat",append=TRUE)

 # End of keywords file
 linea <- paste("HKLOUT ",mtzout)
 cat(linea,file="pointless_keywords.dat",append=TRUE)
 
 # Run program
 stringa <- "pointless < pointless_keywords.dat"
 exepointless <- system(stringa,intern=TRUE,ignore.stderr=TRUE)
 

 # Delete keywords file
 if (file.exists("pointless_keywords.dat")) emptyc <- file.remove("pointless_keywords.dat")

 # Clean directory from unnecessary files
 for (a in files_to_rm) if (file.exists(a)) file.remove(a)

 return(exepointless)
}

#
# Run AIMLESS on specified datasets selection and collect statistical information in a data frame.
merge_datasets <- function(mtz_names,selection,suffix,pointless_keys,aimless_keys,resomin=NULL,resomax=NULL,nref=1)
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

 # Aimless keywords file

 # Default ones
 linea <- "TITLE Scale multiple datasets\n"
 cat(linea,file="aimless_keywords.dat")

 # Check RESOLUTION keyword is already assigned
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
   cnd <- sum(cnd)
  }
 }
 if (length(aimless_keys) == 0) cnd <- 0
 if (cnd == 0)
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

 # Reference file for indexing
 hklref <- indata[nref,1]

 # Run POINTLESS and AIMLESS
 mtz_list <- indata[selection,1]
 sele <- indata[selection,3]
 cat("Merging multiple mtz into a single mtz ...\n")
 exemerge <- merge_mtzs(mtz_list=mtz_list,selection=sele,mtzout="merged.mtz",pointless_keys=pointless_keys,hklref=hklref)
 fatal_error <- grep("#-------------",exemerge,fixed=TRUE)
 if (length(fatal_error) == 0 &
     length(grep("FATAL ERROR message:",exemerge,fixed=TRUE)) == 0 &
     length(exemerge) != 0)
 {
  # Rename merged.mtz
  linea <- paste(suffix[1],paste("merged_",suffix[2],".mtz",sep=""),sep="/")
  file.copy(from="merged.mtz",to=linea)

  # Rename POINTLESS log
  log_file <- paste(suffix[1],paste("pointless_",suffix[2],".log",sep=""),sep="/")
  for (linea in exemerge)
  {
   linea <- paste(linea,"\n",sep="")
   cat(linea,file=log_file,append=TRUE)
  }

  # Prepare command line for AIMLESS
  cat("Running AIMLESS on the merged file ...\n")
  stringa <- sprintf("aimless HKLIN merged.mtz HKLOUT sTrAnO.mtz < aimless_keywords.dat")
  exeaimless <- system(stringa,intern=TRUE, ignore.stderr = TRUE)

  # Run AIMLESS
  nexeaimless <- length(exeaimless)
  if (length(grep("End of aimless job,", exeaimless[(nexeaimless - 1)], fixed = TRUE)) == 0)
  {
   Mstats <- data.frame(Rmeas=c(NA,NA,NA),Rpim=c(NA,NA,NA),Completeness=c(NA,NA,NA),Multiplicity=c(NA,NA,NA),
                        LowRes=c(NA,NA,NA),HighRes=c(NA,NA,NA))

   # AIMLESS log
   log_file <- paste(suffix[1],paste("aimless_",suffix[2],".log",sep=""),sep="/")
   for (linea in exeaimless)
   {
    linea <- paste(linea,"\n",sep="")
    cat(linea,file=log_file,append=TRUE)
   }
 
   # Clean directory from unnecessary files
   dircontents <- system("ls",intern=TRUE)
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

   # Remove files produced to merge mtz's
   dircontents <- system("ls",intern=TRUE)
   idx <- grep("rebatch",dircontents,fixed=TRUE)
   for (i in idx) file.remove(dircontents[i])
   idx <- grep("mtzdump",dircontents,fixed=TRUE)
   for (i in idx) file.remove(dircontents[i])
   idx <- grep("pointless",dircontents,fixed=TRUE)
   for (i in idx) file.remove(dircontents[i])
   #idx <- grep("aimless_keywords",dircontents,fixed=TRUE)
   #for (i in idx) file.remove(dircontents[i])
   idx <- grep("reference",dircontents,fixed=TRUE)
   for (i in idx) file.remove(dircontents[i])
   idx <- grep("merged.mtz",dircontents,fixed=TRUE)
   file.remove(dircontents[idx])
 
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

  # Rename sTrAnO.mtz write AIMLESS log file
  linea <- paste(suffix[1],paste("scaled_",suffix[2],".mtz",sep=""),sep="/")
  file.copy(from="sTrAnO.mtz",to=linea)
  log_file <- paste(suffix[1],paste("aimless_",suffix[2],".log",sep=""),sep="/")
  for (linea in exeaimless)
  {
   linea <- paste(linea,"\n",sep="")
   cat(linea,file=log_file,append=TRUE)
  }

  # Remove files produced to merge mtz's
  dircontents <- system("ls",intern=TRUE)
  idx <- grep("rebatch",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (i in idx) file.remove(dircontents[i])
  idx <- grep("mtzdump",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (i in idx) file.remove(dircontents[i])
  idx <- grep("pointless",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (i in idx) file.remove(dircontents[i])
  idx <- grep("reference",dircontents,fixed=TRUE)
  if (length(idx) != 0) for (i in idx) file.remove(dircontents[i])
  idx <- grep("merged.mtz",dircontents,fixed=TRUE)
  if (length(idx) != 0) file.remove(dircontents[idx])
 
  # Remove files produced to scale mtz's
  dircontents <- system("ls",intern=TRUE)
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
  idx <- grep("sTrAnO.mtz",dircontents,fixed=TRUE)
  if (length(idx) != 0) file.remove(dircontents[idx])

  return(list(Mstats,exeaimless))
 }
 else
 {
  Mstats <- data.frame(Rmeas=c(NA,NA,NA),Rpim=c(NA,NA,NA),Completeness=c(NA,NA,NA),Multiplicity=c(NA,NA,NA),
                       LowRes=c(NA,NA,NA),HighRes=c(NA,NA,NA))
  
  return(list(Mstats,exemerge))
 }
}




######################################################################################################################
######################################################################################################################

# Main


# Load content of previous R run
load("BLEND.RData",.GlobalEnv)

# To avoid warning messages set warn to a negative value 
options(warn = -1)

# Retrieve value from command line
args <- commandArgs(trailingOnly=TRUE)
cat("Files to be combined:  ")
cat(args)
cat("\n")
cat("\n")
cln <- as.numeric(args)

# If one or more datasets have been discarded during the analysis mode, stop
tmp <- match(cln,maindf$cn)
if (sum(is.na(tmp)) > 0) stop("One or more input datasets have been discarded by BLEND during the analysis run, or have never been processed.")

# Find out whether the series of datasets in this exact order has been already processed
gidx <- c()                  # In case there's no GROUPS.info file
outdir <- "./combined_files"
if (file.exists(outdir))
{
 dircontents <- system(paste("ls ",outdir,sep=""),intern=TRUE)

 # Check whether number of files is OK in combined_files directory
 ltmp <- sapply(dircontents,FUN=strsplit,split=".",fixed=TRUE)
 tmp <- unlist(ltmp)
 lidx <- which(tmp == "log")
 midx <- which(tmp == "mtz")
 if (length(midx) != length(lidx))
 {
  messaggio <- "Number of POINTLESS and AIMLESS logs in combined_files directory does not correspond to number of reflection files."
  warning(messaggio)
 }

 # Extract serial number of groups already existing in file GROUPS.info
 mtmp <- paste(outdir,"/GROUPS.info",sep="")
 if (!file.exists(mtmp)) warning("File GROUPS.info is missing in directory combined_files\n")
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
  if (length(lidx) != 2*length(glist)) warning("Number of groups described in file GROUPS.info does not match number of files in directory combined_files")
  for (i in 1:length(glist))
  {
   if (length(glist[[i]]) == length(cln))
   {
    somma <- sum(glist[[i]] == cln)
    if (somma == length(cln)) 
    {
     messaggio <- "You have already merged and scaled this group:  "
     mtmp <- ""
     for (nn in cln) mtmp <- paste(mtmp,nn," ",sep="")
     messaggio <- paste(messaggio,mtmp,"\n",sep="")
     stop(messaggio)
    }
   }
  }
 }

 # Load content of MERGING_STATISTICS.info, if it exists
 mtmp <- paste(outdir,"/MERGING_STATISTICS.info",sep="")
 if (!file.exists(mtmp)) warning("File MERGING_STATISTICS.info is missing in directory combined_files\n")
 if (file.exists(mtmp))
 {
  merging_statistics_info <- scan(mtmp,what="character",sep="\n",quiet=TRUE)
  idx <- grep("Rmeas",merging_statistics_info,fixed=TRUE)
  nn <- length(merging_statistics_info)-idx
  if (nn != length(gidx)) warnings("MERGING_STATISTICS.info has information contrasting with GROUPS_info")
 }
}
if (!file.exists(outdir)) dir.create(outdir)

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
  idxref <- as.integer(tmp[[1]][jdx[length(jdx)]])
 }
 messaggio <- paste("Reference dataset used in case alternate indexing is needed: ",idxref,"\n",sep="")
 cat(messaggio)
 cat("\n")

 # AIMLESS
 idxs <- idxe
 idxe <- length(contents)
 if ((idxe-idxs) > 0) for (line in contents[(idxs+1):idxe]) aimless_keys <- c(aimless_keys,line)
}

# Produce scaled dataset
ftail <- sprintf("%03d",(length(gidx)+1))
suffix <- c(outdir,ftail)
tmp <- merge_datasets("FINAL_list_of_files.dat",selection=cln,suffix,pointless_keys,aimless_keys,
                      resomax=resosel,nref=idxref)
cat(" Statistics for this group:\n")

# Change row names for display purpose
rownames(tmp[[1]]) <- c("Overall","InnerShell","OuterShell")
print(tmp[[1]])
if (length(tmp[[2]]) == 0) warning("No results could be produced for this group due to a problem with POINTLESS")
if (length(tmp[[2]]) != 0 & is.na(tmp[[1]][1,1])) warning("No result could be produced for this group due to a problem with AIMLESS")

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
  linea <- sprintf("   Group Number      Rmeas       Rpim    Completeness    Multiplicity   Lowest Resolution     Highest Resolution\n")
  cat(linea,file=mtmp,append=T)
  linea <- sprintf("            %3d   %8.5f   %8.5f          %6.2f          %6.2f             %7.3f                %7.3f\n",
                   (length(gidx)+1),tmp[[1]][1,1],tmp[[1]][1,2],tmp[[1]][1,3],tmp[[1]][1,4],tmp[[1]][1,5],tmp[[1]][1,6])
  cat(linea,file=mtmp,append=TRUE)
 }
 if (exists("merging_statistics_info"))
 {
  mtmp <- paste(outdir,"/MERGING_STATISTICS.info",sep="")
  linea <- sprintf("            %3d   %8.5f   %8.5f          %6.2f          %6.2f             %7.3f                %7.3f\n",
                   (length(gidx)+1),tmp[[1]][1,1],tmp[[1]][1,2],tmp[[1]][1,3],tmp[[1]][1,4],tmp[[1]][1,5],tmp[[1]][1,6])
  cat(linea,file=mtmp,append=TRUE)
 }
}

# Final lines
#cat("\n")
#cat("##################################################################\n")
#cat("##################################################################\n")
#cat("##################################################################\n")
#cat("  The following files have been created in directory combined_files:\n")
#cat("GROUPS.info                           : ascii file including details (file path, serial number and batch used) for each input group;\n")
#cat("MERGING_STATISTICS.info               : ascii file including a table with merging statistics for each merged and/or scaled file;\n")
#cat("merged_001.mtz, merged_002.mtz, etc   : mtz files prior to scaling. One for each group;\n")
#cat("scaled_001.mtz, scaled_002.mtz, etc   : mtz scaled files. One for each group;\n")
#cat("aimless_001.log, aimless_002.log, etc : log files, one for each AIMLESS job.\n")
#cat("\n")

# Exit without saving
q(save = "no", status = 0, runLast = FALSE)
