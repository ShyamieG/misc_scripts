# AUTHOR: SSG
#
# This script accepts a PLINK-formatted binary file set and outputs an equivalent file set with all duplicate variants (i.e. having the same genomic position) omitted. This script keeps only the first of any duplicated variants it finds.
# This script requires 1) A binary file set name (no extension)
#
# USAGE:
# $ Rscript output_nonDups_bim.R Pagani_2015

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("One argument must be supplied", call.=FALSE)
} else if (length(args)>1) {
  stop("Too many arguments", call.=FALSE)
}

filename <- args[1]

options(stringsAsFactors=F)

cat("reading .bim file...\n")
bim <- read.table(paste(filename, ".bim", sep=""))
new.bim <- cbind(bim, paste(bim$V1, bim$V4, sep=":")); colnames(new.bim)[ncol(new.bim)] <- "pos"

cat("replacing duplicate ID names...\n")
dot.IDs <- which(new.bim$V2==".")
new.bim[dot.IDs, 2] <- new.bim[dot.IDs, "pos"]
i=1; while (any(duplicated(new.bim$V2))) {
    dup.IDs <- new.bim[which(duplicated(new.bim$V2)),2]
    first.match <- match(dup.IDs, new.bim$V2)
    new.bim[first.match, 2] <- paste(new.bim[first.match, 2],i,sep=".")
    i=i+1
}

cat("overwriting .bim file...\n")
write.table(new.bim[,-ncol(new.bim)], paste(filename, ".bim", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
bim <- new.bim

cat("finding and writing out non-duplicate variants...\n")
bim <- bim[match(unique(bim$pos), bim$pos),]
write.table(bim[,-ncol(bim)], paste("nonDupVars_", filename, sep=""), quote=F, row.names=F, col.names=F, sep="\t")

cat("removing duplicates from binary file set...\n")
CMD <- paste("plink2 --allow-no-sex --silent --memory 4000 --threads 1 --bfile ", filename, " --extract nonDupVars_", filename, " --make-bed --out ", filename, "_dupsRemoved", sep="")
system(CMD)