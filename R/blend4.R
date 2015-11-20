#***********************************************************************************************************
#***********************************************************************************************************
#********* blend4.R
#*********
#********* Copyright (C) 2014 Diamond Light Source & Imperial College London
#*********
#********* Authors: James Foadi & Gwyndaf Evans
#*********
#********* This code is distributed under the BSD license, a copy of which is
#********* included in the root directory of this package.
#************************************************************************************************************
#************************************************************************************************************
# Fourth R module of program BLEND (later to be replaced by C++ code)
# It relies on files produced by blend.cpp, blend1.R and blend2.R.
#
#


#
# Functions used in this module
#

# 
# Main function to build the annotated dendrogram
#
annotate_tree <- function(clN,N,Tree,mStats,groups1,macropar,file_name,aLCVs)
{
 # Modify tree distances for better display
 dtmp <- simplify_heights(Tree,groups1)

 # Generate list of clusters
 lista <- find_Nclusters_down(clN,N,dtmp)
 lidx <- lista[[1]]
 lidx2 <- lista[[2]]

 # Change list into vector without NA's
 ulidx <- na.omit(unlist(lidx))

 # Find new nodes coordinates
 nodesxy <- find_nodes_coords(dtmp,groups1,macropar$cn)

 # Data frame with nodes coordinates
 Dnodesxy <- as.data.frame(nodesxy)[ulidx,]

 # Extent of dendrogram details along X
 clns <- groups1[[clN]]
 tmpidx <- match(as.character(clns),dtmp$labels)
 idx <- which(dtmp$order %in% tmpidx)
 xlim <- c(idx[1],idx[length(idx)])

 # y (up and down) extension of clusters window
 DyU <- 0.5
 DyD <- 0.5

 # xlim amd ylim for clusters window
 ylim <- c(min(Dnodesxy$y)-DyD,max(Dnodesxy$y)+DyU)

 # Turn hclust type dendrogram into modern dendrogram object
 ndend <- as.dendrogram(dtmp)

 # Print cluster as zoomed in specified window (PNG)
 fname <- paste(file_name,"png",sep="")
 png(file=fname,width=800,height=1000)
 par(mar=c(5,2,3,4))
 plot(ndend,xlim=xlim,ylim=ylim,yaxt='n',leaflab="none",edgePar=list(lwd=3))
 if (!is.na(aLCVs[1])) 
 {
  points(Dnodesxy,pch=15,cex=6.7,col="grey")
  text(Dnodesxy$x,Dnodesxy$y,labels=ulidx,cex=1.5,col="white",adj=c(0.5,-0.08))
 }
 if (is.na(aLCVs[1])) 
 {
  points(Dnodesxy,pch=16,cex=6,col="grey")
  text(Dnodesxy$x,Dnodesxy$y,labels=ulidx,cex=1.5,col="white")
 }
 labels <- c()
 rownum <- as.integer(rownames(mStats))
 uulidx <- match(ulidx,rownum)
 for (rmeas in mStats$Rmeas[uulidx])
 {
  if (!is.na(rmeas)) tmp <- sprintf("%.3f",rmeas)
  if (is.na(rmeas)) tmp <- ""
  labels <- c(labels,tmp)
 }
 text(Dnodesxy$x,Dnodesxy$y,labels=labels,cex=1,col=4,pos=1,offset=2.0)
 labels <- c()
 for (cmpl in mStats$Completeness[uulidx])
 {
  if (!is.na(cmpl)) tmp <- sprintf("%5.1f%%",cmpl)
  if (is.na(cmpl)) tmp <- ""
  labels <- c(labels,tmp)
 }
 text(Dnodesxy$x,Dnodesxy$y,labels=labels,cex=1,col=3,adj=c(1.3,-2.5))
 labels <- c()
 for (cc12 in mStats$CC12[uulidx])
 {
  if (!is.na(cc12)) 
  {
   tmp <- sprintf("%.1f",cc12)
   tmp <- paste(tmp,"\uc5")
  }
  if (is.na(cc12)) tmp <- ""
  labels <- c(labels,tmp)
 }
 text(Dnodesxy$x,Dnodesxy$y,labels=labels,cex=1,col=2,adj=c(-0.3,-2.5))
 greyL <- rgb(0.5,0.5,0.5)
 points(idx,rep(ylim[1],times=length(idx)),type="h",col=greyL)

 # Print aLCV_values only for "a" and "aDO" modes (where we pass NA, rather than aLCV_values)
 if (!is.na(aLCVs[1]))
 {
  labels <- c()
  for (ii in uulidx)
  {
   tmp <- sprintf("%7.2f",aLCVs[ii])
   labels <- c(labels,tmp)
  }
  text(Dnodesxy$x,Dnodesxy$y,labels=labels,cex=0.9,col=4,adj=c(0.6,1.65))
 }

 # Print data set number only if tree branch ends in data set, not in cluster
 edges <- find_nodes_edges(dtmp,groups1,macropar)
 text <- c()
 labels <- c()
 for (ii in 1:length(ulidx))
 {
  uu <- ulidx[ii]
  ee <- edges[[uu]]
  if (ee[1]%%1 < 0.000001) 
  {
   eles <- -dtmp$merge[uu,1]
   if (eles == dtmp$order[as.integer(ee[1])])
   {
    text <- c(text,ee[1])
    labels <- c(labels,dtmp$labels[eles])
   }
  }
  if (ee[2]%%1 < 0.000001) 
  {
   eles <- -dtmp$merge[uu,2]
   if (eles == dtmp$order[as.integer(ee[2])])
   {
    text <- c(text,ee[2])
    labels <- c(labels,dtmp$labels[eles])
   }
  }
 }
 if (length(text) > 0) mtext(text=labels,side=1,at=text,col=2,line=1,cex=1.0)
 dev.off()

 # Print cluster as zoomed in specified window (PS)
 fname <- paste(file_name,"ps",sep="")
 postscript(file=fname,paper = "a4",horizontal=FALSE)
 par(mar=c(5,2,3,4))
 plot(ndend,xlim=xlim,ylim=ylim,yaxt='n',leaflab="none",edgePar=list(lwd=3))
 if (!is.na(aLCVs[1])) 
 {
  points(Dnodesxy,pch=15,cex=6.7,col="grey")
  text(Dnodesxy$x,Dnodesxy$y,labels=ulidx,cex=1.5,col="white",adj=c(0.5,-0.08))
 }
 if (is.na(aLCVs[1])) 
 {
  points(Dnodesxy,pch=16,cex=8,col="grey")
  text(Dnodesxy$x,Dnodesxy$y,labels=ulidx,cex=2.3,col="white")
 }
 labels <- c()
 rownum <- as.integer(rownames(mStats))
 uulidx <- match(ulidx,rownum)
 for (rmeas in mStats$Rmeas[uulidx])
 {
  if (!is.na(rmeas)) tmp <- sprintf("%.3f",rmeas)
  if (is.na(rmeas)) tmp <- ""
  labels <- c(labels,tmp)
 }
 text(Dnodesxy$x,Dnodesxy$y,labels=labels,cex=0.7,col=4,pos=1,offset=1.7)
 labels <- c()
 for (cmpl in mStats$Completeness[uulidx])
 {
  if (!is.na(cmpl)) tmp <- sprintf("%5.1f%%",cmpl)
  if (is.na(cmpl)) tmp <- ""
  labels <- c(labels,tmp)
 }
 text(Dnodesxy$x,Dnodesxy$y,labels=labels,cex=0.7,col=3,adj=c(1.2,-3.9))
 labels <- c()
 for (cc12 in mStats$CC12[uulidx])
 {
  if (!is.na(cc12)) 
  {
   tmp <- sprintf("%.1f",cc12)
   tmp <- paste(tmp,"\uc5")
  }
  if (is.na(cc12)) tmp <- ""
  labels <- c(labels,tmp)
 }
 text(Dnodesxy$x,Dnodesxy$y,labels=labels,cex=0.7,col=2,adj=c(-0.2,-3.9))
 greyL <- rgb(0.5,0.5,0.5)
 points(idx,rep(ylim[1],times=length(idx)),type="h",col=greyL)

 # Print aLCV_values only for "a" and "aDO" modes (where we pass NA, rather than aLCV_values)
 if (!is.na(aLCVs[1]))
 {
  labels <- c()
  for (ii in uulidx)
  {
   tmp <- sprintf("%7.2f",aLCVs[ii])
   labels <- c(labels,tmp)
  }
  text(Dnodesxy$x,Dnodesxy$y,labels=labels,cex=0.9,col=4,adj=c(0.6,1.65))
 }

 # Print data set number only if tree branch ends in data set, not in cluster
 edges <- find_nodes_edges(dtmp,groups1,macropar)
 text <- c()
 labels <- c()
 for (ii in 1:length(ulidx))
 {
  uu <- ulidx[ii]
  ee <- edges[[uu]]
  if (ee[1]%%1 < 0.000001) 
  {
   eles <- -dtmp$merge[uu,1]
   if (eles == dtmp$order[as.integer(ee[1])])
   {
    text <- c(text,ee[1])
    labels <- c(labels,dtmp$labels[eles])
   }
  }
  if (ee[2]%%1 < 0.000001) 
  {
   eles <- -dtmp$merge[uu,2]
   if (eles == dtmp$order[as.integer(ee[2])])
   {
    text <- c(text,ee[2])
    labels <- c(labels,dtmp$labels[eles])
   }
  }
 }
 if (length(text) > 0) mtext(text=labels,side=1,at=text,col=2,line=1,cex=1.0)
 dev.off()

 return(list(ulidx,Dnodesxy))
}

#
# Function to change heights of dendrogram nodes so that clusters with same number of
# elements are at the same height. Heights proceeds in steps of 1
#
simplify_heights <- function(Tree,groups1)
{
 # Create vector with number of elements in each cluster
 ltmp <- lapply(groups1,length)
 ulist <- unlist(ltmp)

 # Assign order number to increasingly complex clusters
 idss <- unique(ulist)
 Tree$height <- match(ulist,idss)

 return(Tree)
}

#
# Find elements and cluster under a node, up to N levels downward
#
find_Nclusters_down <- function(clN,N,Tree)
{
 # Initialise list
 lidx <- list(clN)

 # Fill list
 for (i in 2:N)
 {
  tmp <- c()
  for (j in 1:length(lidx[[i-1]]))
  {
   if (is.na(lidx[[i-1]][j]))
   {
    tmp <- c(tmp,NA,NA)
   }
   if (!is.na(lidx[[i-1]][j]))
   {
    tmp <- c(tmp,Tree$merge[lidx[[i-1]][j],])
   }
  }
  tmp[tmp < 0] <- NA
  lidx[[i]] <- tmp
 }

 # A second list without NA's might be useful
 lidx2 <- list(clN)
 for (i in 2:length(lidx))
 {
  tmp <- c()
  for (j in 1:length(lidx2[[i-1]]))
  {
   if (lidx2[[i-1]][j] > 0) tmp <- c(tmp,Tree$merge[lidx2[[i-1]][j],])
   if (lidx2[[i-1]][j] < 0) tmp <- c(tmp,lidx2[[i-1]][j])
  }
  lidx2[[i]] <- tmp
 }

 # Prune at level when all elements are NA
 newlidx <- list()
 for (i in 1:N)
 {
  if (sum(lidx[[i]],na.rm=TRUE) > 0) newlidx[[i]] <- lidx[[i]]
 }

 return(list(newlidx=newlidx,lidx2=lidx2))
}

#
# Find coordinates of tree edges at all nodes
#
find_nodes_edges <- function(uTree,groups1,macropar)
{
 # Number of clusters
 nclusters <- length(uTree$merge[,1])

 # Unique heights
 levels <- sort(unique(uTree$height))

 # Find coordinates of cluster nodes
 nodesxy <- find_nodes_coords(uTree,groups1,macropar$cn)

 # Cycle through unique heights
 LList <- list()
 RList <- list()
 for (i in 1:nclusters)
 {
  if (uTree$height[i] == 1)
  {
   xcentre <- nodesxy$x[i]
   ltmp <- as.integer(xcentre-0.5)
   rtmp <- as.integer(xcentre+0.5)
   LList[[i]] <- ltmp
   RList[[i]] <- rtmp
  }
  if (uTree$height[i] > 1)
  {
   eles <- uTree$merge[i,]
   if (eles[1] < 0)
   {
    idx <- which(uTree$order == -eles[1])
    LList[[i]] <- idx
   }
   if (eles[1] > 0) LList[[i]] <- nodesxy$x[eles[1]]
   if (eles[2] < 0)
   {
    idx <- which(uTree$order == -eles[2])
    RList[[i]] <- idx
   }
   if (eles[2] > 0) RList[[i]] <- nodesxy$x[eles[2]]
  }
 }
 edges <- list(L=unlist(LList),R=unlist(RList))

 # Better a list of couples
 tmp <- list()
 for (i in 1:length(edges$L))
 {
  tmp[[i]] <- c(edges$L[i],edges$R[i])
 }
 edges <- tmp

 return(edges)
}


######################################################################################################################
######################################################################################################################

# Main

# Find out if running on Windows
ostuff <- Sys.info()
rwin = FALSE
if (ostuff["sysname"] == "Windows") rwin = TRUE

# To avoid warning messages set warn to a negative value 
options(warn = -1)

# Retrieve value from command line
args <- commandArgs(trailingOnly=TRUE)
cat("\n")
cat(paste("Input string:   ",paste(args,collapse=" ")))
n <- length(args)
cat("\n")
cat("\n")

# Decide what type of graphics job it is required, based on first argument.
#
# D   :   Build PNG annotated dendrograms
#

# Exit if none of the graphics modes is met
if (args[1] != "D" & args[1] != "DO")
{
 msg <- sprintf("Type %s in command line is unrecognised\n",args[1])
 cat(msg)
 q(save = "no", status = 1, runLast = FALSE)
}

# D: Build PNG annotated dendrograms
if (args[1] == "D" | args[1] == "DO")
{

 # Load contents of previous R run from blend1.R and blend2.R
 if (file.exists("BLEND0.RData")) load("BLEND0.RData",.GlobalEnv)
 if (file.exists("BLEND.RData")) load("BLEND.RData",.GlobalEnv)
 if (!file.exists("BLEND.RData") & !file.exists("BLEND0.RData"))
 {
  msg <- "Neither file BLEND.RData nor BLEND0.RData found.\n"
  cat(msg)
  q(save = "no", status = 1, runLast = FALSE)
 }

 # Change groups according to where information has been loaded from
 if (is.list(groups[[1]]))
 {
  groups1 <- groups[[1]]
  groups2 <- groups[[2]]
 }
 if (!is.list(groups[[1]]))
 {
  groups1 <- groups
  groups2 <- NULL
 }

 # Check if synthesis mode has been executed
 if (file.exists("BLEND.RMergingStatistics")) 
 {
  load("BLEND.RMergingStatistics",.GlobalEnv)

  # We also need to pass NA, rather than aLCV_values, as these won't be displayed
  aLCVs <- NA
 }
 if (!file.exists("BLEND.RMergingStatistics") | args[1] == "DO")
 {
  # Create a mergingStatistics data frame with all NA's, because synthesis has not been carried out yet

  # Number of dendrogram nodes
  ncls <- length(groups1)

  # Vectors including low resolutions and high resolutions of clusters
  if (!is.null(groups2))
  {
   vecLow <- c()
   vecHigh <- c()
   for (i in 1:ncls) 
   {
    vecLow <- c(vecLow,groups2[[i]][1])
    vecHigh <- c(vecHigh,groups2[[i]][2])
   }
  }
  if (is.null(groups2))
  {
   vecLow <- rep(NA,times=ncls)
   vecHigh <- rep(NA,times=ncls)
  }

  # Data frame mergingStatistics
  mergingStatistics <- data.frame(Rmeas=rep(NA,times=ncls),Rpim=rep(NA,times=ncls),
                                  Completeness=rep(NA,times=ncls),                  # Later to be changed according to new code
                                  Multiplicity=rep(NA,times=ncls),                  # Later to be changed according to new code
                                  LowRes=vecLow,HighRes=vecHigh,
                                  CC12=rep(NA,times=ncls),                          # Later to be changed according to new code
                                  Mn2=rep(NA,times=ncls))                           # Later to be changed according to new code

  # We also need to pass aLCV_values, as these will be displayed
  aLCVs <- aLCV_values
 }

 # Default cluster is top cluster and number of clusters levels downward is 2
 if (n == 1) 
 {
  clN <- length(groups1) 
  msg <- "WARNING! Cluster number unassigned. Default to top cluster.\n"
  cat(msg)
  cat("\n")
  N <- 2
  msg <- "WARNING! Dendrogram level unassigned. Default to level 2.\n"
  cat(msg)
  cat("\n")
 }
 if (n == 2)
 {
  tmp <- as.integer(args[2])
  if (!is.na(tmp)) 
  {
   clN <- tmp
   if (clN > length(groups1))
   {
    msg <- "WARNING! Cluster number higher than top cluster. Default to top cluster.\n"
    cat(msg)
    cat("\n")
    clN <- length(groups1)
   }
   if (clN < 1)
   {
    msg <- "WARNING! Cluster number less than one. Default to top cluster.\n"
    cat(msg)
    cat("\n")
    clN <- length(groups1)
   }
  }
  if (is.na(tmp)) 
  {
   clN <- length(groups1)
   msg <- "WARNING! Unrecognised cluster. Default to top cluster.\n"
   cat(msg)
   cat("\n")
  }
  N <- 2
  msg <- "WARNING! Dendrogram level unassigned. Default to level 2.\n"
  cat(msg)
  cat("\n")
 }
 if (n >= 3)
 {
  tmp <- as.integer(args[2])
  if (!is.na(tmp)) 
  {
   clN <- tmp
   if (clN > length(groups1))
   {
    msg <- "WARNING! Cluster number higher than top cluster. Default to top cluster.\n"
    cat(msg)
    cat("\n")
    clN <- length(groups1)
   }
   if (clN < 1)
   {
    msg <- "WARNING! Cluster number less than one. Default to top cluster.\n"
    cat(msg)
    cat("\n")
    clN <- length(groups1)
   }
  }
  if (is.na(tmp)) 
  {
   clN <- length(groups1)
   msg <- "WARNING! Unrecognised cluster. Default to top cluster.\n"
   cat(msg)
   cat("\n")
  }
  tmp <- as.integer(args[3])
  if (!is.na(tmp)) 
  {
   N <- tmp
   if (N < 2)
   {
    msg <- "WARNING! Dendrogram level smaller than 2. Default to level 2.\n"
    cat(msg)
    cat("\n")
    N <- 2
   }
  }
  if (is.na(tmp)) 
  {
   N <- 2
   msg <- "WARNING! Unrecognised dendrogram level Default to level 2.\n"
   cat(msg)
   cat("\n")
  }
 }

 # Create directory "graphics", if it does not exists
 if (!file.exists("./graphics")) dir.create("./graphics")

 # Create annotated dendrogram PNG file base name
 if (!is.na(aLCVs[1]))
 {
  file_name <- sprintf("./graphics/aLCV_annotated_dendrogram_cluster_%03d_level_%03d.",clN,N)
 }
 if (is.na(aLCVs[1]))
 {
  file_name <- sprintf("./graphics/stats_annotated_dendrogram_cluster_%03d_level_%03d.",clN,N)
 }

 # Create annotated dendrograms (png and ps)
 tmp <- annotate_tree(clN,N,npar.hc_ward,mergingStatistics,groups1,macropar,file_name,aLCVs)

 # Normal termination
 q(save = "no", status = 0, runLast = FALSE)
}
