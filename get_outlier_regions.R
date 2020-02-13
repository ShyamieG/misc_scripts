# AUTHOR: SSG
# This script accepts a list of outlier SNPs (3 column file: chromosome, SNP position, and count)
# It returns a list of regions
# The required argument is: Input (the name of the file)
# The optional arguments are: MaxRange - the minimum physical distance (in bp) required between two SNPs to start a new region - default is 3 Mbp
#                           : Buffer - the physical distance to add to the beginning and end of a region - default is 0.25 Mbp
# USAGE:
# $ Rscript get_outlier_regions.R GERMLINE2_IBDdepth_IBDoutliers.txt 3e6 2.5e5

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Defaults
maxRange <- 3e6
buffer <- 2e5

# Error handling and updating default values
if (length(args)<1) {
  stop("At least one argument must be supplied", call.=FALSE)
} else if (length(args)>=2) {
  maxRange <- as.numeric(args[2])
} else if (length(args)==3) {
  buffer <- as.numeric(args[3])
} else if (length(args)>3) {
  stop("Too many arguments", call.=FALSE)
}

FILE <- read.table(args[1], stringsAsFactors=F)

find.IBDoutlier.range <- function(pileup.outliers, maxRange, buffer) {
  outlier.regions <- as.data.frame(matrix(ncol=4, nrow=0))
    colnames(outlier.regions) <- c("CHR", "START", "END", "TYPE")
  # Split into chromosomes and sort the outlier values (going along the chromosome)
  chrs <- unique(pileup.outliers[,1])
  outliers.by.chr <- list()
  for (i in 1:length(chrs)) {
    outliers.by.chr[[i]] <- pileup.outliers[pileup.outliers[,1]==chrs[i],]
    outliers.by.chr[[i]] <- outliers.by.chr[[i]][order(outliers.by.chr[[i]][,2]),2]
  }
  for (i in 1:length(chrs)) {
    # If there is only one outlier SNP, define as the region +/- the buffer around that SNP
    if (length(outliers.by.chr[[i]])==1) {
      outlier.range <- c(outliers.by.chr[[i]][1] - buffer, outliers.by.chr[[i]][1] + buffer)
      outlier.regions[nrow(outlier.regions)+1,] <- c(chrs[i], outlier.range, "depth_outlier")
    }
    # If there is more than one outlier SNP, initiate a new region that begins with the first SNP
    if (length(outliers.by.chr[[i]])>1) {
      outlier.range <- c()
      outlier.range[1] <- outliers.by.chr[[i]][1] - buffer
      # Find endpoints of the region
      endpoint.inds <- which(outliers.by.chr[[i]][1:((length(outliers.by.chr[[i]]))-1)] - outliers.by.chr[[i]][-1] < -maxRange)
      endpoint.inds <- c(endpoint.inds, length(outliers.by.chr[[i]]))
      for (n in 1:length(endpoint.inds)) {
        if (n==length(endpoint.inds)) {
          outlier.range[2] <- outliers.by.chr[[i]][endpoint.inds[n]] + buffer
          outlier.regions[nrow(outlier.regions)+1,] <- c(chrs[i], outlier.range, "depth_outlier")
        } else {
          outlier.range[2] <- outliers.by.chr[[i]][endpoint.inds[n]] + buffer
          outlier.regions[nrow(outlier.regions)+1,] <- c(chrs[i], outlier.range, "depth_outlier")
          outlier.range[1] <- outliers.by.chr[[i]][endpoint.inds[n]+1] - buffer
        }
      }
    }
  }
  return(outlier.regions)
}

results <- find.IBDoutlier.range(FILE, maxRange, buffer)
write.table(results, "", quote=F, col.names=F, row.names=F, sep=" ")