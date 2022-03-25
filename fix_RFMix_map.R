# AUTHOR: SSG
# This script accepts a genomic map file containing 3 columns - chromosome, base pair position, and centiMorgan position (for RFMix)
# This script returns a genomic map file with runs of decreasing positions removed (this is necessary because liftOver sometimes results in these runs)
# The required arguments are: Input map file - full path to the genomic map file to be fixed
#			      Offset value - number of rows in the map scan to look before and after a run to optimize the window to remove in order to fix the map
#			      Output map file - full path to the desired output file
# USAGE:
# $ Rscript fix_RFMix_map.R input_genome.map 100 output_genome.map

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("At least three arguments must be supplied", call.=FALSE)
} else if (length(args)>3) {
  stop("Too many arguments", call.=FALSE)
}

options(stringsAsFactors=F)
map <- read.table(args[1]);colnames(map) <- c("chr", "bp", "cM")
offset <- as.numeric(args[2])
outpath <- args[3]

find.opt.win <- function(dat, run, metric, offset) {
  check.pos <- function(pos, dat, metric) {
    start <- as.numeric(pos[1])
    end <- as.numeric(pos[2])
    if (dat[end,metric] > dat[start,metric]) {
      return(t(as.data.frame(c(start, end, end-start))))
    }
  }
  starts <- max(1, (min(run)-offset)):(min(run)-1)
  ends <- (max(run)+1):min((max(run)+offset), nrow(dat))
  combos <- expand.grid(starts, ends)
  if (dat[max(ends),metric] < dat[min(starts),metric]) {
    # if end of chromosome is still too low, just delete it
    opt.result <- c(min(run), max(ends), max(ends)-min(run));names(opt.result) <- c("start","end", "diff")
  } else {
    # otherwise, find optimal window to remove
    results <- as.data.frame(matrix(unlist(apply(X=combos, MARGIN=1, FUN=check.pos, dat=dat, metric=metric)), ncol=3, byrow=T))
    colnames(results) <- c("start","end","diff")
    opt.result <- unlist(results[which.min(results$diff),])
  }
  return(opt.result)
}

for (chr in 1:22) {
  dat <- map[map$chr==paste0("chr", chr),]
  # note how many rows were in original data
  old.nrow<- nrow(dat)
  # if chromosome map doesn't start at lowest genomic distance, remove all rows up to the lowest genomic distance
  for (metric in c("bp", "cM")) {
    if (dat[1,metric] != min(dat[,metric])) {
      dat <- dat[-(1:(which.min(dat[,metric])-1)),]
      cat(paste("Removed a decreasing stretch of", length((1:(which.min(dat[,metric])-1))), metric, "values at the start of chromosome", chr, "\n"))
    }
  }
  for (metric in c("bp", "cM")) {
    diff.points <- which(diff(dat[,metric]) < 0)
    if (length(diff.points) > 0) {
      decr.runs <- split(seq_along(diff.points), cumsum(c(0, diff(diff.points) > 1)))
    } else {
      decr.runs <- NULL
    }
    if (length(decr.runs) > 1) {
      cat(paste0("On chromosome ", chr, ", found ", length(decr.runs), " stretches of decreasing ", metric, " values "))
    } else if (length(decr.runs) == 1) {
      cat(paste0("On chromosome ", chr, ", found ", length(decr.runs), " stretch of decreasing ", metric, " values "))
    }
    if (length(diff.points) > 1) {
      cat(paste0("spanning ", length(diff.points), " rows\n"))
    } else if (length(diff.points) == 1) {
      cat(paste0("spanning ", length(diff.points), " row\n"))
    }
    while (length(decr.runs) > 0) {
      run <- diff.points[decr.runs[[1]]]
      # remove runs that are at the start or end of chromosome
      if (run[1] == 1 | run[length(run)] >= nrow(dat)) {
        dat <- dat[-run,]
      } else {
        opt.win.to.remove <- find.opt.win(dat, run, metric, offset)
        dat <- dat[-((opt.win.to.remove["start"]+1):(opt.win.to.remove["end"]-1)),]
      }
      diff.points <- which(diff(dat[,metric]) < 0)
      if (length(diff.points) > 0) {
        decr.runs <- split(seq_along(diff.points), cumsum(c(0, diff(diff.points) > 1)))
      } else {
        decr.runs <- NULL
      }
    }
  }
  new.nrow <- nrow(dat)
  if (chr == 1) {
    new.map <- dat
  } else {
    new.map[(nrow(new.map)+1):(nrow(new.map)+nrow(dat)),] <- dat
  }
  cat(paste("Finished chromosome", chr, "-", old.nrow-new.nrow, "rows removed from map\n"))
}

cat("Writing out new map...")
write.table(new.map, outpath, quote=F, row.names=F, col.names=F, sep="\t")
cat(" done!\n")