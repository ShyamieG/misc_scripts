# AUTHOR: SSG
#
# Calculate population by population Weir-Cockerham Fst (adapted from Laura R Botigue's script)
# This script requires 1) The file name for a set of binary plink files
#
# It also accepts an optional code file that tells the script with individuals belong to which population (three columns, the first having the population code, and the next two being the first two columns of the .fam file)
# If this is not provided, it is assumed that the first column of the .fam file contains the population codes
#
# This script returns a matrix of population by population Weir-Cockerham Fst values
#
# USAGE:
# $ Rscript calc_pop_WC_Fst.R [PLINK file name] [code file]]
#
# ADDITIONAL NOTES:
# The .fam file MUST contain only -9 in the phenotype column

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("At least one argument must be supplied", call.=FALSE)
} else if (length(args)>2) {
  stop("Too many arguments", call.=FALSE)
}

FAM <- read.table(paste(args[1], "fam", sep="."), stringsAsFactors=F)
  SKIP=F
  if (sum(FAM[,6]) != nrow(FAM)*-9) {
    warning("Phenotype column of .fam contains non '-9' values")
    SKIP=T
  }

cat("Calculating global SNP frequencies for the dataset\n")
system2("plink2", paste("--bfile", args[1], "--silent --freq --out", args[1]))
FRQ <- read.table(paste(args[1], "frq", sep="."), stringsAsFactors=F, header=T)
MAJ.ALLELES <- FRQ[,3]

if (length(args)==1) {
  POPS <- unique(FAM$V1)
  codefile <- FAM[,c(1,1:6)]
} else if (length(args)==2) {
  codefile <- read.table(args[2], stringsAsFactors=F)
  # Put in same order as .fam file
  codefile <- merge(codefile, FAM, by.x=2:3, by.y=1:2)[,c(3,1,2,4:7)]
  if (nrow(codefile)!=nrow(FAM)) {
    stop("Code file must contain 1 row for every individual in the .fam file", call.=FALSE)
  }
  POPS <- unique(codefile$V1)
}

Fst.MATRIX <- as.data.frame(matrix(nrow=length(POPS), ncol=length(POPS)))
  colnames(Fst.MATRIX) <- POPS; rownames(Fst.MATRIX) <- colnames(Fst.MATRIX)
# Number of individuals per population
Fst_N <- c()
# Hets per SNP per population
HETS.By.POP <- as.data.frame(matrix(nrow=nrow(FRQ), ncol=length(POPS)))
  colnames(HETS.By.POP) <- POPS
  rownames(HETS.By.POP) <- FRQ[,"SNP"]
  
# Freqs per SNP per population
FRQS.By.POP <- as.data.frame(matrix(nrow=nrow(FRQ), ncol=length(POPS)))
  colnames(FRQS.By.POP) <- POPS
  rownames(FRQS.By.POP) <- FRQ[,"SNP"]

# Loop over every population to calculate SNP stats
cat("Calculating population level stats per SNP\n")
for (P in POPS) {
  KEEP <- codefile[grep(P, codefile[,1]), 2:7]
  Fst_N <- c(Fst_N, nrow(KEEP))
  write.table(KEEP, "temp_keep.txt", quote=FALSE, col.names=FALSE, row.names=FALSE, sep=" ")
  
  # Calculate number of heterozygotes per SNP within the population
  system2("plink2", paste("--bfile", args[1], "--silent --keep temp_keep.txt --hardy --out temp"))
  # Extract het count per SNP
  HWE.file <- read.table("temp.hwe", stringsAsFactors=F, header=TRUE)
  hets <- unlist(lapply(strsplit(as.character(HWE.file[,6]), "/"), "[", 2))
  if (SKIP==F) {
    HETS.By.POP[,P] <- hets
  } else if (SKIP==T) {
    HETS.By.POP[,P] <- hets[seq(1, length(hets), 3)]
  }
  
  # Calculate allele frequency at each SNP within the population
  system2("plink2", paste("--bfile", args[1], "--silent --keep temp_keep.txt --freq --out temp"))
  # Extract frq per SNP
  FRQ.file <- read.table("temp.frq", stringsAsFactors=F, header=TRUE)
  frqs <- as.numeric(FRQ.file[,5])
  # For loci where the major allele has flipped, use the inverse frequency
  allele.flip <- which(FRQ.file[,3] != MAJ.ALLELES)
  frqs[allele.flip] <- 1 - as.numeric(FRQ.file[allele.flip,5])
  
  FRQS.By.POP[,P] <- frqs
  
  cat(paste("  ", P, "done\n", sep=" "))
}
system2("rm", "temp.hwe temp.frq temp_keep.txt temp.nosex temp.log")

cat("Writing tables\n")
write.table(HETS.By.POP, "Hets_by_Pop.txt", quote=FALSE, col.names=FALSE, row.names=TRUE, sep=" ")
write.table(FRQS.By.POP, "Freqs_by_Pop.txt", quote=FALSE, col.names=FALSE, row.names=TRUE, sep=" ")

cat("Calculating population pairwise Fst\n")

Fst.MATRIX[POPS[length(POPS)], POPS[length(POPS)]] <- 0
for (i in 1:(length(POPS)-1) ) {
  Fst.MATRIX[POPS[i], POPS[i]] <- 0
  for (j in (i+1):length(POPS)) {
    cat(paste("  Starting ", POPS[i], " by ", POPS[j], "\n", sep=""))
    
    n = Fst_N[c(i, j)] # no. of inds in each pop being considered
    n_avg = mean(n)
    r=2 ##Number of populations
    
    sum_nc <- 0
    for (k in 1:r) {
      sum_nc = sum_nc + ((n[k]^2) /(r*n_avg)) # sample size squared / 2 * mean pop size 
    }
    nc = (r*n_avg - sum_nc) / (r-1) # (pop # x mean - val from above) / (pop # - 1)
    
    # Find # of hets and allele frequencies per pop
    Fst.hets <- HETS.By.POP[,c(POPS[i],POPS[j])]
    Fst.freqs <- FRQS.By.POP[,c(POPS[i],POPS[j])]
    NA.SNPs <- unique(c(which(is.na(Fst.freqs[,1])), which(is.na(Fst.freqs[,2]))))
    
    # Remove NA values
    if ( length(NA.SNPs) != 0 ) {
      Fst.hets <- Fst.hets[-NA.SNPs,]
      Fst.freqs <- Fst.freqs[-NA.SNPs,]
    }
    
    p_avg <- 0
    for (k in 1:r) {
      p_avg = p_avg + (n[k] * as.numeric(Fst.freqs[,k]) / (r*n_avg)) # sum of sample size times allele frequency of a population divided by the number of populations times the average population size
    }
    
    h_avg <- 0
    s2 <- 0
    for (k in 1:r) {
      s2 = s2 + (n[k] * (as.numeric(Fst.freqs[, k]) -p_avg)^2 /((r-1)*n_avg))
      h_avg = h_avg + (as.numeric(Fst.hets[,k]) / (r*n_avg)) 
    }
    
    # Calculate values a, b, and c    
    a = (n_avg/nc) * (s2 - ((1/(n_avg-1)) * (p_avg*(1-p_avg) - ((r-1) / r)*s2 - 0.25*h_avg)))
	b = (n_avg / (n_avg-1)) * (p_avg*(1-p_avg) - ((r-1) / r)*s2 - (((2*n_avg - 1) / (4*n_avg)) * h_avg))
	c = 0.5 * h_avg
    
    b[b<0] <- 0
    
    Fst.MATRIX[POPS[i], POPS[j]] <- sum(a)/(sum(a,b,c))
    Fst.MATRIX[POPS[j], POPS[i]] <- Fst.MATRIX[POPS[i], POPS[j]]
  }
}

cat("Writing Fst matrix\n")
write.table(Fst.MATRIX, "Fst_matrix.txt", quote=F, sep=" ")