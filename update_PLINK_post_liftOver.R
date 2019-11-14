# AUTHOR: SSG
#
# This script updates PLINK files after liftOver has been conducted on the corresponding UCSC format .bed file
# This script requires 1) The old PLINK .bim file, generated from the old PLINK .bim file
#                      2) The new (UCSC format) .bed file, generated from liftOver
#                      3) A list of unmapped SNPs that are missing in the new .bed file (modified from liftOver output, should contain only list of missing SNPs, chrom and chromEnd)
#
# USAGE:
# $ Rscript update_PLINK_post_liftOver.R [old bim] [new bed] [missing SNP list]

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("Four arguments must be supplied", call.=FALSE)
} else if (length(args)>3) {
  stop("Too many arguments", call.=FALSE)
}

options(stringsAsFactors=F, header=F)
cat("reading files...\n")
oldbim <- read.table(args[1])
newbed <- read.table(args[2])
missingSNPs <- read.table(args[3])

# Extract PLINK binary file name
bim.filename <- unlist(strsplit(args[1], split=".", fixed=T))
bim.filename <- paste(bim.filename[-match("bim", bim.filename)], collapse=".")

cat("checking files...\n") # some checks for inconsistencies
if ((nrow(newbed) + nrow(missingSNPs)) != nrow(oldbim)) {
	stop("Total number of variants in unlifted file plus the new .bed file does not equal the number of variants in the original bed file")
}

old.var.counts <- table(paste(oldbim$V1, oldbim$V4))
dup.vars <- names(old.var.counts)[which(old.var.counts > 1)]

if (length(dup.vars) > 1) {
	stop("Duplicate variants exist in .bim file - remove these with output_nonDups_bim.R before proceeding")
}

cat("processing...\n")
# Remove 'chr' from chromosome name, if it is there
missingSNPs$V1 <- unlist(strsplit(missingSNPs$V1, split="chr"))[seq(2, nrow(missingSNPs)*2, 2)]

# Create the new .bim file
newbed.2 <- newbed[,c(1,3)]
newbed.2$V1 <- unlist(strsplit(newbed$V1, split="chr"))[seq(2, nrow(newbed)*2, 2)]

newbim <- oldbim[-match(paste(missingSNPs$V1, missingSNPs$V3), paste(oldbim$V1, oldbim$V4)),]
newbim[,c(1,4)] <- newbed.2

cat("generating new files...\n")
write.table(newbim, "SNPs_to_keep_from_old_bim", row.names=F, col.names=F, quote=F)
CMD <- paste("plink2 --bfile ", bim.filename, " --extract SNPs_to_keep_from_old_bim --silent --memory 4000 --threads 1 --make-bed --out ", paste(bim.filename, "liftedOver", sep="_"), sep="")
system(CMD)
system("rm SNPs_to_keep_from_old_bim")
write.table(newbim, paste(bim.filename, "_liftedOver.bim", sep=""), row.names=F, col.names=F, quote=F)
