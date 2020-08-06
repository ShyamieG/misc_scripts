dropin <- read.table("dropInSamplesInfo", stringsAsFactors=F, header=T)
PCAs <- dropin$sample_ID
n_eigs=10

for (i in PCAs) {
    eigs.file <- read.table(paste(i, "_PCA.eigs", sep=""), row.names=1)[,-(n_eigs+1)]
    rownames(eigs.file) <- unlist(strsplit(rownames(eigs.file), split=":"))[seq(2, nrow(eigs.file)*2, 2)]
    write.table(eigs.file, paste(i, "_PCA.csv", sep=""), quote=F, row.names=T, col.names=F, sep=",")
}