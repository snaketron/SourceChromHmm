

# Description:
# Computes HDI given a vector, taken "Doing Bayesian Analysis"
getHdi <- function(vec, hdi.level) {
  sortedPts <- sort(vec)
  ciIdxInc <- floor(hdi.level * length(sortedPts))
  nCIs = length(sortedPts) - ciIdxInc
  ciWidth = rep(0 , nCIs)
  for (i in 1:nCIs) {
    ciWidth[i] = sortedPts[i + ciIdxInc] - sortedPts[i]
  }
  HDImin = sortedPts[which.min(ciWidth)]
  HDImax = sortedPts[which.min(ciWidth) + ciIdxInc]
  HDIlim = c(HDImin, HDImax)
  return(HDIlim)
}


# Description:
# Get the names of the 22 mouse chromosomes
getChrs <- function() {
  return(c(paste("chr", 1:19, sep = ''), "chrX", "chrY", "chrMt"))
}


# Description:
# Get the uset ChIP marks
getMarks <- function() {
  marks <- c("H3K4me3", "H3K9ac", "H3K27ac", 
             "H3K36me3", "H3K9me3", "H3K27me3")
  return(marks)
}


# Description:
# Get all genes
getGenes <- function(dir) {
  genes <- read.csv(file = dir, sep = "\t", as.is = T, header = T)
  genes$tss <- ifelse(test = genes$Strand == 1, yes = genes$Gene.start, no = genes$Gene.end)
  genes$chr <- paste("chr", genes$Chr, sep = '')
  genes$chr[genes$chr == "chrMT"] <- "chrMt"
  genes$key <- paste(genes$Gene.name, genes$Gene.stable.ID, sep = '_')
  genes$bp.start <- genes$Gene.start
  genes$bp.stop <- genes$Gene.end
  genes$bp.mid <- genes$Gene.start+(genes$Gene.end-genes$Gene.start)/2
  return(genes)
}



getMm10File <- function(dir) {
  out <- read.csv(file = dir, header = F, sep = "\t", as.is = T)
  colnames(out) <- c("chr", "bp.start", "bp.stop")
  out$tss <- out$bp.start
  return (out)
}
