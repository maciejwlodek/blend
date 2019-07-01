#***********************************************************************************************************
#***********************************************************************************************************
#********* blend0.R
#*********
#********* Copyright (C) 2014 Diamond Light Source & Imperial College London
#*********
#********* Authors: James Foadi & Gwyndaf Evans
#*********
#********* This code is distributed under the BSD license, a copy of which is
#********* included in the root directory of this package.
#************************************************************************************************************
#************************************************************************************************************
# Zeroeth R module of program BLEND (later to be replaced by C++ code)
# It relies on files produced by blend.cpp.
#
#


# Load libraries and auxiliary files
require(MASS, quietly = TRUE, warn.conflicts = FALSE)  ## For rlm (robust regression)

# Functions

printf <- function(...) cat(sprintf(...))

distance <- function(p1, p2) {
    diff <- p1-p2
    sum_sq <- 0
    for(i in 1:length(diff)) {
        sum_sq <- sum_sq + diff[i]^2
    }
    return(sqrt(sum_sq))
}
print_clusters <- function(datasets, k, labels) {
    cat("CLUSTERS.txt\n")
    for(i in 1:k) {
        printf("Cluster %d\n", i)
        for(j in 1:nrow(datasets)) {
            if(labels[j] == i) {
                cat(j)
                cat(" ")
            }
        }
        cat("\n")
    }
}

update_centroid <- function(datasets, index, labels, centroid) {
    dim <- ncol(datasets)
    sum <- 0
    for(i in 1:(dim-1)) sum <- c(sum, 0)
    counter <- 0
    for(i in 1:length(labels)) {
        if(labels[i] == index) {
            sum <- sum + datasets[i, ]
            counter <- counter + 1
        }
    }
    if(counter == 0) return(centroid)
    centroid <- sum / counter
    return(centroid)
}
update_labels <- function(datasets, centroids, labels, k) {
    for(x in 1:nrow(datasets)) {
        dat <- datasets[x, ]
        min_dist <- .Machine$double.xmax
        for(i in 1:k) {
            cent <- centroids[i, ]
            d <- distance(dat, cent)
            if(d < min_dist) {
                min_dist <- d
                labels[x] <- i
            }
        }
    }
    return(labels)
}
kmeans_clust <- function(datasets, k, num_iters) {
    npoints <- nrow(datasets)
    dim <- ncol(datasets)
    rand_indices <- sample.int(npoints, k)
    centroids <- c()
    for(i in 1:k) {
        centroids <- c(centroids, datasets[rand_indices[i], 2:7])
    }
    centroids <- matrix(as.vector(centroids), nrow=k, ncol=dim, byrow=TRUE)
    labels <- vector("list", npoints)
    labels <- update_labels(datasets, centroids, labels, k)
    j <- 0
    converged <- FALSE
    while(j < num_iters && converged == FALSE) {
        converged <- TRUE
        for(i in 1:k) {
            tmp <- centroids[i, ]
            centroids[i, ] <- update_centroid(datasets, i, labels, tmp)
            if(length(setdiff(tmp, centroids[i, ])) != 0 || length(setdiff(centroids[i, ], tmp)) != 0) {
                converged <- FALSE
            }
        }
        labels <- update_labels(datasets, centroids, labels, k)
        j <- j+1
    }
    if(converged == FALSE) {
        printf("did not converge in %d iterations!\n", num_iters)
    }
    return(labels)
}
listGroups <- function(datasets, k, labels) {
    groups <- list()
    for(i in 1:k) {
        group_i <- list()
        for(j in 1:nrow(datasets)) {
            if(labels[j] == i) {
                group_i <- c(group_i, j)
            }
        }
        group_i <- as.integer(group_i)
        if(length(group_i) > 0) groups <- c(groups, list(group_i))
    }
#    print(groups)
    return(groups)
}

# To extract a list with datasets to be merged, as many as the nodes of a dendrogram (htree)
treeToList2 <- function(htree)
{
 # Function to list all groups corresponding to each node in a dendrogram.
 # "htree" is the cluster object (from hclust).


 groups <- list()
 itree <- dim(htree$merge)[1]
 for (i in 1:itree)
 {
  il <- htree$merge[i,1]
  ir <- htree$merge[i,2]
  if (il < 0) lab1 <- htree$labels[-il]
  if (ir < 0) lab2 <- htree$labels[-ir]
  if (il > 0) lab1 <- groups[[il]]
  if (ir > 0) lab2 <- groups[[ir]]
  lab <- c(lab1,lab2)
  lab <- as.integer(lab)
  groups <- c(groups,list(lab))
 }

 return(groups)
}

#
# Reduce data dimensionality and standardize for cell parameters
nparCell <- function(data,cn,cumprob)
{
 # Normalize data
 model <- prcomp(data,scale=TRUE)
 smod <- summary(model)
 # Choose minimal number of variables that give enough statistical variability
 if (length(model$x[1,]) == 1) npar <- model$x
 if (length(model$x[1,]) > 1)
 {
  idx <-  which(smod$importance[3,] > cumprob)
  npar <- model$x[,1:idx[1]]
 }
 npar <- as.matrix(npar)
 rownames(npar) <- cn

 return(npar)
}

# Functions to calculate cell variation in angstroms
.dcl <- function(a,b,c,aa,bb,gg)
{
# Input: cell parameters. Output: values useful to all crystallographic
# calculations. These are:
# 1) sa = sin(alpha), sb = sin(beta), sc = sin(gamma)
# 2) ca = cos(alpha), cb = cos(beta). cc = cos(gamma)
# 3) sides of reciprocal cell: ar = a*, br = b*, cr = c*
# 4) sines of angles of reciprocal cell: sar = sin(alpha*), sbr = sin(beta*), scr = sin(gamma*)
# 5) cosines of angles of reciprocal cell: car = cos(alpha*), cbr = cos(beta*), ccr = cos(gamma*)
# 6) Volume of unit cell: V
 aa <- aa*pi/180
 bb <- bb*pi/180
 gg <- gg*pi/180
 sa <- sin(aa)
 sb <- sin(bb)
 sc <- sin(gg)
 ca <- cos(aa)
 cb <- cos(bb)
 cc <- cos(gg)

 # To avoid NaN generated by rounding off errors, use for cell-derived quantities formulas
 # derived previously by computationa crystallographers
 sang <- 0.5*(aa+bb+gg)
 V2 <- sqrt(sin(sang-aa)*sin(sang-bb)*sin(sang-gg)*sin(sang))
 V <- 2*a*b*c*V2
 ar <- b*c*sa/V
 br <- a*c*sb/V
 cr <- a*b*sc/V
 car <- (cb*cc-ca)/(sb*sc)
 cbr <- (ca*cc-cb)/(sa*sc)
 ccr <- (ca*cb-cc)/(sa*sb)
 sar <- sqrt(1-car*car)
 sbr <- sqrt(1-cbr*cbr)
 scr <- sqrt(1-ccr*ccr)
 l <- c(sa,sb,sc,ca,cb,cc,ar,br,cr,sar,sbr,scr,car,cbr,ccr,V)
 names(l) <- c("SIN_ALPHA","SIN_BETA","SIN_GAMMA","COS_ALPHA","COS_BETA","COS_GAMMA","A*","B*","C*","SIN_ALPHA*","SIN_BETA*","SIN_GAMMA*",
               "COS_ALPHA*","COS_BETA*","COS_GAMMA*","V")

 return(l)
}

# Given cell parameters, return the inverse of cell orthogonalization matrix (First choice in Giacovazzo's book)
.triclinic_to_orthogonal_01 <- function(a,b,c,aa,bb,cc)
{
 lp <- .dcl(a,b,c,aa,bb,cc)
 m_1 <- matrix(c(a,0,0,b*lp[6],b*lp[3],0,c*lp[5],-c*lp[2]*lp[13],1/lp[9]),nrow=3,ncol=3,byrow=TRUE)

 return(m_1)
}

# Given cell parameters, return the inverse of cell orthogonalization matrix (MOSFLM choice, second choice in Giacovazzo's book)
.triclinic_to_orthogonal_02 <- function(a,b,c,aa,bb,cc)
{
 lp <- .dcl(a,b,c,aa,bb,cc)
 m_1 <- matrix(c(1/lp[7],-lp[15]/(lp[7]*lp[12]),a*lp[5],0,1/(lp[8]*lp[12]),b*lp[4],0,0,c),nrow=3,ncol=3,byrow=TRUE)

 return(m_1)
}

# Compute matrix of cross-cells max distances
evaluateMaxChange <- function(maindf)
{
 n <- length(maindf[,2])
 dMat <- matrix(nrow=n,ncol=n)
 for (i in 1:n)
 {
  for (j in 1:n)
  {
   M1 <- .triclinic_to_orthogonal_01(maindf[i,2],maindf[i,3],maindf[i,4],maindf[i,5],maindf[i,6],maindf[i,7])
   M2 <- .triclinic_to_orthogonal_01(maindf[j,2],maindf[j,3],maindf[j,4],maindf[j,5],maindf[j,6],maindf[j,7])
   v1 <- M1%*%matrix(c(1,1,1),ncol=1)
   v2 <- M2%*%matrix(c(1,1,1),ncol=1)
   dd <- sqrt((v1[1,1]-v2[1,1])^2+(v1[2,1]-v2[2,1])^2+(v1[3,1]-v2[3,1])^2)
   dMat[i,j] <- dd
  }
 }

 return(dMat)
}

# Functions to evaluate max change in linear dimensions for cell parameters

# Diagonals along 3 independent faces of unit cell
faceDiagonals <- function(a,b,c,aa,bb,cc)
{
 dab <- sqrt(a^2+b^2-2*a*b*cos(pi-cc*pi/180))
 dac <- sqrt(a^2+c^2-2*a*c*cos(pi-bb*pi/180))
 dbc <- sqrt(b^2+c^2-2*b*c*cos(pi-aa*pi/180))

 return(c(dab,dac,dbc))
}

distRatio <- function(v)
{
 # Length of vector
 n <- length(v)

 # Create n X n matrix
 m <- matrix(ncol=n,nrow=n)

 # Double loop to fill the matrix
 for (i in 1:n)
 {
  for (j in 1:n)
  {
   dd <- abs(v[i]-v[j])
   mv <- min(v[i],v[j])
   m[i,j] <- (dd/mv)*100
  }
 }

 return(m)
}

adistRatio <- function(v)
{
 # Length of vector
 n <- length(v)

 # Create n X n matrix
 m <- matrix(ncol=n,nrow=n)

 # Double loop to fill the matrix
 for (i in 1:n)
 {
  for (j in 1:n)
  {
   m[i,j] <- abs(v[i]-v[j])
  }
 }

 return(m)
}

# Get matrix indices
get_indices <- function(N,n)
{
 # N is the serial number and n the size of n X n matrix
 i <- N%%n
 if (i == 0) i <- n
 j <- N%/%n+1
 if (i == n) j <- j-1
 if (j > n) j <- n

 return(c(i,j))
}

maxRatio <- function(macropar,idx)
{
 # Cell parameters
 cpar <- macropar[idx,2:7]

 # Number of datasets
 n <- length(cpar[,1])
 
 # 3 diagonal lengths for all datasets
 dab <- c()
 dac <- c()
 dbc <- c()
 for (i in 1:n)
 {
  tmp <- faceDiagonals(cpar[i,1],cpar[i,2],cpar[i,3],cpar[i,4],cpar[i,5],cpar[i,6])
  dab <- c(dab,tmp[1])
  dac <- c(dac,tmp[2])
  dbc <- c(dbc,tmp[3])
 }

 # Calculate maxRatio matrix for the 3 diagonals vectors and extract max value for each matrix
 #mab <- max(distRatio(dab))
 mab <- max(adistRatio(dab))
 #iab <- which(distRatio(dab) == max(distRatio(dab)))
 iab <- which(adistRatio(dab) == max(adistRatio(dab)))
 ijab <- get_indices(iab[1],n)
 #sab <- adistRatio(dab)[iab]
 sab <- distRatio(dab)[iab]
 #mac <- max(distRatio(dac))
 mac <- max(adistRatio(dac))
 #iac <- which(distRatio(dac) == max(distRatio(dac)))
 iac <- which(adistRatio(dac) == max(adistRatio(dac)))
 ijac <- get_indices(iac[1],n)
 #sac <- adistRatio(dac)[iac]
 sac <- distRatio(dac)[iac]
 #mbc <- max(distRatio(dbc))
 mbc <- max(adistRatio(dbc))
 #ibc <- which(distRatio(dbc) == max(distRatio(dbc)))
 ibc <- which(adistRatio(dbc) == max(adistRatio(dbc)))
 ijbc <- get_indices(ibc[1],n)
 #sbc <- adistRatio(dbc)[ibc]
 sbc <- distRatio(dbc)[ibc]
 vv <- c(sab[1],sac[1],sbc[1])
 ij <- matrix(c(ijab,ijac,ijbc),ncol=3)
 imall <- which(c(mab,mac,mbc) == max(c(mab,mac,mbc)))[1]
 Mpar <- vv[imall]
 ijs <- ij[,imall]
 cns <- macropar[idx[ijs],"cn"]

 #return(c(max(mab,mac,mbc),Mpar,cns[1],cns[2]))
 return(c(Mpar,max(mab,mac,mbc),cns[1],cns[2]))
}


find_nodes_coords <- function(clst, clns, cn)
{
 # Number of objects
 nobj <- length(clst$merge[,1]) + 1

 # x coordinates for individual objects
 xobj <- clns[[nobj - 1]]

 # Go through nodes formation
 x <- numeric(nobj - 1)
 y <- numeric(nobj - 1)
 for (i in 1:length(clst$merge[,1]))
 {
  ele1 <- clst$merge[i,1]
  ele2 <- clst$merge[i,2]
  if (ele1 < 0 & ele2 < 0)
  {
   #idx1 <- match(abs(ele1), xobj)
   #idx2 <- match(abs(ele2), xobj)
   idx1 <- match(cn[abs(ele1)], xobj)
   idx2 <- match(cn[abs(ele2)], xobj)
   x[i] <- 0.5 * (idx1 + idx2)
   y[i] <- clst$height[i]
  }
  if (ele1 < 0 & ele2 > 0)
  {
   #idx1 <- match(abs(ele1), xobj)
   idx1 <- match(cn[abs(ele1)], xobj)
   idx2 <- x[ele2]
   x[i] <- 0.5 * (idx1 + idx2)
   y[i] <- clst$height[i]
  }
  if (ele1 > 0 & ele2 < 0)
  {
   idx1 <- x[ele1]
   #idx2 <- match(abs(ele2), xobj)
   idx2 <- match(cn[abs(ele2)], xobj)
   x[i] <- 0.5 * (idx1 + idx2)
   y[i] <- clst$height[i]
  }
  if (ele1 > 0 & ele2 > 0)
  {
   idx1 <- x[ele1]
   idx2 <- x[ele2]
   x[i] <- 0.5 * (idx1 + idx2)
   y[i] <- clst$height[i]
  }
 }

 return(list(x=x,y=y))
}
####################################################################
####################################################################
######################### Main program #############################
####################################################################
####################################################################

# To avoid warning messages set warn to a negative value 
options(warn = -1)

# Fix global parameters here (default values)
nbin <- 20                    # Number of bins for dynamic and overall Wilson plots
fdrop <- 0.75                 # Intensity decay fraction for datasets affected by radiation damage (keyword RADFRAC)
isigi <- 1.5                  # Signal-to-noise ratio wanted to estimate highest resolution cut
cparweight <- 1.0             # Weight to be given to cluster analysis with cell parameters as descriptors.
                              # Wilson plots descriptors get automatically weight 1-cparweight
k <- 3

#######################################################################################################
#######################################################################################################
######################################## KEYWORDS FROM FILE ###########################################
#######################################################################################################
#######################################################################################################
# Check whether there's a ready keyword file and use it if found
if (file.exists("BLEND_KEYWORDS.dat"))
{
 contents <- scan("BLEND_KEYWORDS.dat", what="character",sep="\n",quiet=TRUE)

 # Find section "BLEND KEYWORDS"
 idxs <- grep("BLEND KEYWORDS",contents,fixed=TRUE)
 idxe <- grep("POINTLESS KEYWORDS",contents,fixed=TRUE)

 # Only extract values if something is included in "BLEND KEYWORDS" section
 if ((idxe-idxs) > 1)
 {
  idxs <- idxs+1
  idxe <- idxe-1

  # Turn values into data frame
  slm <- strsplit(contents[idxs:idxe]," ")
  tmp <- data.frame()
  for (slt in slm)
  {
   jdx <- which(nchar(slt) != 0)
   #if (length(jdx) != 2) stop("Wrongly formatted BLEND_KEYWORDS.dat file")
   tmp <- rbind(tmp,data.frame(I(slt[jdx[1]]),I(slt[jdx[2]])))             # "I()" to avoid data frame turning characters into factors
  }

  for (i in 1:length(tmp[,1]))
  {
   if (as.character(tmp[i,1]) == "NBIN" | as.character(tmp[i,1]) == "nbin") nbin <- as.numeric(tmp[i,2])
   if (as.character(tmp[i,1]) == "RADFRAC" | as.character(tmp[i,1]) == "radfrac") fdrop <- as.numeric(tmp[i,2])
   if (as.character(tmp[i,1]) == "ISIGI" | as.character(tmp[i,1]) == "isigi") isigi <- as.numeric(tmp[i,2])
   if (as.character(tmp[i,1]) == "CPARWT" | as.character(tmp[i,1]) == "cparwt") cparweight <- as.numeric(tmp[i,2])
   if (as.character(tmp[i,1]) == "K" | as.character(tmp[i,1]) == "k") k <- as.numeric(tmp[i,2])
  }

  # Turn indices back to their initial value for following section
  idxs <- idxs-1
  idxe <- idxe+1
 }

 # Not all keywords can be included in existing "BLEND_KEYWORDS.dat". Complete as needed.
 cat("BLEND KEYWORDS\n",file="BLEND_KEYWORDS.dat")
 if ((idxe-idxs) > 1)
 {
  idxs <- idxs+1
  idxe <- idxe-1
  for (i in idxs:idxe) cat(paste(contents[i],"\n",sep=""),file="BLEND_KEYWORDS.dat",append=TRUE)
  nomi <- as.character(tmp[,1])
  if (!("NBIN" %in% nomi) & !("nbin" %in% nomi))
  {
   linea <- sprintf("NBIN      %d\n",as.integer(nbin))
   cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
  }
  if (!("RADFRAC" %in% nomi) & !("radfrac" %in% nomi))
  {
   linea <- sprintf("RADFRAC   %5.3f\n",as.numeric(fdrop))
   cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
  }
  if (!("ISIGI" %in% nomi) & !("isigi" %in% nomi))
  {
   linea <- sprintf("ISIGI     %5.3f\n",as.numeric(isigi))
   cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
  }
  if (!("CPARWT" %in% nomi) & !("cparwt" %in% nomi))
  {
   linea <- sprintf("CPARWT    %5.3f\n",as.numeric(cparweight))
   cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
  }
  if (!("K" %in% nomi) & !("k" %in% nomi)) {
   linea <- sprintf("K         %d\n", as.integer(k))
   cat(linea, file="BLEND_KEYWORDS.dat", append=TRUE)
  }

  # Turn indices back to their initial value for following section
  idxs <- idxs-1
  idxe <- idxe+1
 }

 # If nothing was included in BLEND KEYWORDS section add everything
 if ((idxe-idxs) == 1)
 {
  linea <- sprintf("NBIN      %d\n",as.integer(nbin))
  cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
  linea <- sprintf("RADFRAC   %5.3f\n",as.numeric(fdrop))
  cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
  linea <- sprintf("ISIGI     %5.3f\n",as.numeric(isigi))
  cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
  linea <- sprintf("CPARWT    %5.3f\n",as.numeric(cparweight))
  cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
  linea <- sprintf("K         %d\n", as.integer(k))
  cat(linea,file="BLEND_KEYWORDS.dat", append=TRUE)
 }

 # Append rest of keywords
 for (i in idxe:(length(contents)-1)) cat(paste(contents[i],"\n",sep=""),file="BLEND_KEYWORDS.dat",append=TRUE)
 cat(contents[length(contents)],file="BLEND_KEYWORDS.dat",append=TRUE)
}

# If keyword file is not found write one with default words
if (!file.exists("BLEND_KEYWORDS.dat"))
{
 cat("BLEND KEYWORDS\n",file="BLEND_KEYWORDS.dat")
 linea <- sprintf("NBIN      %d\n",as.integer(nbin))
 cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
 linea <- sprintf("RADFRAC   %5.3f\n",as.numeric(fdrop))
 cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
 linea <- sprintf("ISIGI     %5.3f\n",as.numeric(isigi))
 cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
 linea <- sprintf("CPARWT    %5.3f\n",as.numeric(cparweight))
 cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
 linea <- sprintf("K         %d\n", as.integer(k))

 # Complete with keywords sections for POINTLESS and AIMLESS
 cat("POINTLESS KEYWORDS\n",file="BLEND_KEYWORDS.dat",append=TRUE)
 cat("AIMLESS KEYWORDS\n",file="BLEND_KEYWORDS.dat",append=TRUE)
}

#######################################################################################################
#######################################################################################################
###################################### DESCRIPTORS SECTION ############################################
######################################   CELL PARAMETERS   ############################################
#######################################################################################################
#######################################################################################################
# Load data from file "forR_macropar.dat"
macropar <- read.table("forR_macropar.dat")
names(macropar) <- c("cn","a","b","c","alpha","beta","gamma","mosa","ctoddist","wlength")

# Re-arrange macropar with same order as file NEW_list_of_files.dat
macropar <- macropar[order(macropar$cn),]

# Prepare data frame for cluster analysis
cat("Preparing data for cluster analysis ...\n")
maindf <- macropar[,1:7]

# Decide which macroparameters can be used in cluster analysis
cell_columns <- c()
if (sum(is.na(maindf$a)) == 0)
{
 tmp <- sd(maindf$a)
 if (tmp > 0.000001) cell_columns <- c(cell_columns,2)
}
if (sum(is.na(maindf$b)) == 0)
{
 tmp <- sd(maindf$b)
 if (tmp > 0.000001) cell_columns <- c(cell_columns,3)
}
if (sum(is.na(maindf$c)) == 0)
{
 tmp <- sd(maindf$c)
 if (tmp > 0.000001) cell_columns <- c(cell_columns,4)
}
if (sum(is.na(maindf$alpha)) == 0)
{
 tmp <- sd(maindf$alpha)
 if (tmp > 0.000001) cell_columns <- c(cell_columns,5)
}
if (sum(is.na(maindf$beta)) == 0)
{
 tmp <- sd(maindf$beta)
 if (tmp > 0.000001) cell_columns <- c(cell_columns,6)
}
if (sum(is.na(maindf$gamma)) == 0)
{
 tmp <- sd(maindf$gamma)
 if (tmp > 0.000001) cell_columns <- c(cell_columns,7)
}

# Full list of descriptors used
fullc <- cell_columns

#######################################################################################################
#######################################################################################################
########################################## CLUSTER ANALYSIS ###########################################
#######################################################################################################
#######################################################################################################
# Enough information to start cluster analysis (add bit concerning cluster analysis here)
cat("Cluster analysis initiated ...\n")

# Reduce data dimensionality and standardize
if (length(fullc) != 0) 
{
 nparC <- nparCell(maindf[,fullc],maindf$cn,cumprob=0.95)

 # Build distance matrix for cluster analysis (2 methods)
 nparC.dist <- dist(nparC)
}

# Final distance matrix has contributions from both set of descriptors
if (length(fullc) == 0 | cparweight < 1)
{
 msg <- sprintf("CPARWT parameter cannot be different from 1 in dendrogram-only mode.")
 cat(msg)
 q(save = "no", status = 1, runLast = FALSE)
}

# Cluster analysis

#maindftest <- maindf[,3:4]
#maindftest <- data.matrix(maindftest)
#k <- 3
niters <- 30
#labels <- kmeans_clust(maindftest, k, niters)
#print_clusters(maindftest, k, labels)
labels <- kmeans_clust(nparC, k, niters)
print_clusters(nparC, k, labels)
#groups <- listGroups(maindftest, k, labels)
groups <- listGroups(nparC, k, labels)

# Calculate LCV values for all clusters
LCV_values <- c()
aLCV_values <- c()
LCV_couples <- matrix(nrow=length(groups),ncol=2)
icpls <- 0
for (cn in groups)
{
 icpls <- icpls+1
 idx <- match(cn,macropar$cn)
 #tmp <- maxRatio(macropar[idx,2:7])
 tmp <- maxRatio(macropar,idx)
 LCV_values <- c(LCV_values,tmp[1])
 aLCV_values <- c(aLCV_values,tmp[2])
 minV <- min(tmp[3],tmp[4])
 maxV <- max(tmp[3],tmp[4])
 LCV_couples[icpls,1] <- as.integer(minV)
 LCV_couples[icpls,2] <- as.integer(maxV)
}

# Output merging nodes table to an ascii file
nTable <- "./CLUSTERS.txt"
if (file.exists(nTable)) emptyc <- file.remove(nTable)
linea <- "                               \n"
cat(linea,file=nTable)
linea <- " Cluster     Number of           LCV      aLCV   Furthest    Datasets\n"
cat(linea,file=nTable,append=TRUE)
linea <- "  Number      Datasets                           Datasets          ID\n"
cat(linea,file=nTable,append=TRUE)
linea <- "                               \n"
cat(linea,file=nTable,append=TRUE)
for (i in 1:length(groups))
{
 sorted_groups <- sort(groups[[i]])
 linea <- paste(sprintf("     %03d           %3d       %7.2f %9.2f %5d %4d  ",
                i,length(groups[[i]]),LCV_values[i],aLCV_values[i],LCV_couples[i,1],LCV_couples[i,2]),"  ",
                paste(sorted_groups,collapse=" "),"\n",sep="")
 cat(linea,file=nTable,append=TRUE)
}

cat("Cluster analysis completed!\n")

# Remove "forR_macropar.dat" file
if (file.exists("forR_macropar.dat")) emptyc <- file.remove("forR_macropar.dat")

# Write file with original crystal number
filenames <- read.table("./NEW_list_of_files.dat",as.is=c(1))
emptyc <- file.remove("./NEW_list_of_files.dat")
for (ii in 1:length(filenames[,1]))
{
 stmp <- normalizePath(filenames$V1[ii])
 filenames$V1[ii] <- stmp
}
if (file.exists("FINAL_list_of_files.dat")) emptyc <- file.remove("FINAL_list_of_files.dat")
idx <- match(1:length(filenames[,1]),maindf$cn) 
for (i in 1:length(filenames[,1]))
{
 if (!is.na(idx[i]))
 {
  linea <- sprintf("%-130s %5.0f\n",filenames[i,],maindf$cn[idx[i]])
 }
 if (is.na(idx[i]))
 {
  linea <- sprintf("%-130s %5.0ff\n",filenames[i,],NA)
 }
 cat(linea,file="./FINAL_list_of_files.dat",append=T)
}

# Save image for next run of BLEND
save.image(file="BLEND0.RData")

# Exit without saving
q(save = "no", status = 0, runLast = FALSE)
