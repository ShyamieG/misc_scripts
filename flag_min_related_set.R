# AUTHOR: SSG
# This script accepts the output of plink --genome and a relatedness threshold
# This script returns (a shortest) list of individuals that would need to be removed to elminate all close relationships in the dataset
# The required argument is: Input (the stem of the .genome and .fam files)
# The optional argument is: PI_HAT threshold - the minimum value of PI_HAT that defines a related pair (default 0.125)
# USAGE:
# $ Rscript flag_min_related_set.R CHABU 0.1 > min_relateds.list

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("At least one argument must be supplied", call.=FALSE)
} else if (length(args)==1) {
  pihat.max <- 1/8
} else if (length(args)==2) {
  pihat.max <- args[2]
} else if (length(args)>2) {
  stop("Too many arguments", call.=FALSE)
}

KINSHIP.DAT <- read.table(paste(args[1], "genome", sep="."), stringsAsFactors=F, header=T)

min.related.set <- function(kinship.data, pihat.max) {
  relateds <- c()
  relateds.data <- kinship.data[kinship.data[,"PI_HAT"]>pihat.max,]
  while (nrow(relateds.data)>=1) {
    ID <- names(which.max(table(c(relateds.data[,"IID1"], relateds.data[,"IID2"]))))
    relateds.data <- relateds.data[-c(grep(ID, relateds.data[,"IID1"]), grep(ID, relateds.data[,"IID2"])),]
    relateds <- c(relateds, ID)
  }
  return(relateds)
}

min.related.set(KINSHIP.DAT, pihat.max)
