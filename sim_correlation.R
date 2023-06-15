library(vioplot)

# observing a moderate correlation across array and methyl-seq does not necessarily support the finding of true ancestry-specific differences in DNA methylation

# choose sample size
N <- 211
# choose target correlation coefficient
r <- 0.6
# choose AFR allele frequency
AFR_p <- 0.5
# choose mean of group with the allele and without the allele
g1_mean <- 0.7
g3_mean <- 0.2
# choose standard deviation
SD <- 0.8
# choose group colours
group1col <- rgb(30/255, 140/255, 50/255, 0.7)
group2col <- rgb(200/255, 190/255, 70/255, 0.7)
group3col <- rgb(170/255, 50/255, 120/255, 0.7)

##### ----------------------------------------------

# determine allele frequency at probe SNP
p <- AFR_p * 0.8
# determine groups
group1 <- 1:round(p^2*N)
group2 <- (round(p^2*N)+1):(round(p^2*N)+round(2*(1-p)*p*N))
group3 <- (round(p^2*N)+round(2*(1-p)*p*N)+1):N
# determine mean DNA methylation level for heterozygotes
g2_mean <- mean(c(g1_mean, g3_mean))

# make up some array data with a systematic difference between two groups
set.seed(12345)
array_dat <- c(rnorm(length(group1), mean=g1_mean, sd=SD), rnorm(length(group2), mean=g2_mean, sd=SD), rnorm(length(group3), mean=g3_mean, sd=SD))

# make up some methylation data that is unimodal
set.seed(45678)
mseq_dat0 <- rnorm(N, g2_mean, sd=SD)
# adjust values to exhibit the desired correlation with array_dat
mseq_dat = r*array_dat + sqrt(1-r*r) * mseq_dat0

# effect sizes and p-values
DAT <- as.data.frame(cbind(c(rep(1, length(group1)), rep(2, length(group2)), rep(3, length(group3))), array_dat, mseq_dat));colnames(DAT)[1] <- "group"
DAT[,2:3] <- (DAT[,2:3] - min(DAT[,2:3]))/(max(DAT[,2:3])-min(DAT[,2:3])) # re-scale values to range between 0 and 1, like DNA methylation beta values
# groups are significantly different by array
X <- lm(DAT$array_dat ~ DAT$group)
cat(paste0("SNP effect size of ", signif(X$coefficients[[2]], 3), " from simulated array data (p = ", signif(summary(X)$coefficients[2,4], 3), ")\n"))
# groups are NOT significantly different by methylation sequencing
Y <- lm(DAT$mseq_dat ~ DAT$group)
cat(paste0("SNP effect size of ", signif(Y$coefficients[[2]], 3), " from simulated methylation data (p = ", signif(summary(Y)$coefficients[2,4], 3), ")\n"))

# we can achieve the target correlation without reproducing the between-group differences
m <- rbind(c(1, 2),c(3, 3));layout(m) # set up plot format
vioplot(DAT$array_dat[DAT$group==1], DAT$array_dat[DAT$group==2], DAT$array_dat[DAT$group==3], main="simulated array data", col=c(group1col, group2col, group3col), ylab="DNA methylation value", xlab="genotype at probe SNP", ylim=c(0,1))
vioplot(DAT$mseq_dat[DAT$group==1], DAT$mseq_dat[DAT$group==2], DAT$mseq_dat[DAT$group==3], main="simulated methyl-seq data", col=c(group1col, group2col, group3col), ylab="DNA methylation value", xlab="genotype at probe SNP", ylim=c(0,1))
cor <- cor(DAT$array_dat, DAT$mseq_dat)
plot(DAT$mseq_dat[DAT$group==1], DAT$array_dat[DAT$group==1], col=group1col, pch=19, xlim=c(0,1), ylim=c(0,1), xlab="simulated methyl-seq data", ylab="simulated array data")
points(DAT$mseq_dat[DAT$group==2], DAT$array_dat[DAT$group==2], col=group2col, pch=19)
points(DAT$mseq_dat[DAT$group==3], DAT$array_dat[DAT$group==3], col=group3col, pch=19)
text(min(DAT[,2:3]), max(DAT[,2:3])-(max(DAT[,2:3])-min(DAT[,2:3]))*0.04, paste("correlation =", signif(cor,3)), pos=4)
