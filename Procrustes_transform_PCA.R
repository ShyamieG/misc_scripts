## AUTHOR: SSG
## This script performs the Procrustes transformation on the dropped in samples within the Procrustes_transform_PCA.py script
#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
library(vegan)
options(stringsAsFactors=F)

n_eigs <- as.numeric(args[1])
outfile <- "Procrustes_PCA_out.csv"

## Determine 'base' PCA to use - based on which individual had the most SNPs available
SNP_counts <- read.table("dropInSamplesInfo", header=F)
ranked_inds <- SNP_counts[order(SNP_counts[,2], decreasing=T),1]
base_ind = ranked_inds[1]

## Write out eigenvalues for base PCA
Header <- readLines(paste(base_ind, "_PCA.eigs", sep=""), n=1)
Header <- unlist(strsplit(Header, split=" ")); Header <- Header[which(Header!="")]
write.table(t(c("", Header[-1])), file=outfile, quote=F, sep=",", row.names=F, col.names=F)

## Write out base PC values, including best drop in sample
base_PC_vals <- read.table(paste(base_ind, "_PCA.eigs", sep=""), row.names=1)[,-(n_eigs+1)]
rownames(base_PC_vals) <- unlist(strsplit(rownames(base_PC_vals), split=":"))[seq(2, nrow(base_PC_vals)*2, 2)]
write.table(base_PC_vals, file=outfile, quote=F, sep=",", col.names=F)
refOnly_base_PC_vals <- base_PC_vals[-match(base_ind, rownames(base_PC_vals)),]

## Perform Procrustes transformation on the rest of the drop in samples
for (ind in ranked_inds[-1]) {
    cat(paste("Performing Procrustes transform on", ind, "\n"))
    dropin_PC_vals <- read.table(paste(ind, "_PCA.eigs", sep=""), row.names=1)[,-(n_eigs+1)]
    rownames(dropin_PC_vals) <- unlist(strsplit(rownames(dropin_PC_vals), split=":"))[seq(2, nrow(dropin_PC_vals)*2, 2)]
    refOnly_dropin_PC_vals <- dropin_PC_vals[-match(ind, rownames(dropin_PC_vals)),]
    refOnly_dropin_PC_vals <- refOnly_dropin_PC_vals[rownames(refOnly_base_PC_vals),]
    # Determine scale, rotation, translation to minimize sum of squares w/ target points
    Proc_transform=procrustes(X=as.matrix(refOnly_base_PC_vals), Y=as.matrix(refOnly_dropin_PC_vals))
    # adjust PCs for drop in sample by parameters determined in previous step
    output=predict(Proc_transform, dropin_PC_vals[match(ind, rownames(dropin_PC_vals)),])
    write.table(output, file=outfile, append=T, col.names=F, row.names=T, quote=F, sep=",")
}

## fin