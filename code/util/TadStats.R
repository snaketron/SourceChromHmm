source("code/util/GeneralUtilities.R")

# Description:
# Gene is in tad if its TSS is in the TAD
getTssTadStats <- function(tad.dir, genes.data) {
  tad <- read.csv(file = tad.dir, sep = '\t', header = F, as.is = T)
  colnames(tad) <- c("chr", "bp.start", "bp.stop")
  tad$bp.mid <- tad$bp.start+(tad$bp.stop-tad$bp.start)/2
  tad$mb.length <- (tad$bp.stop-tad$bp.start)/10^6
  tad$chr <- factor(x = tad$chr, levels = getChrs())
  tad$range <- paste(tad$bp.start, tad$bp.stop, sep = '-')
  tad$ID <- paste(tad$chr, 1:nrow(tad), tad$range, sep = '.')
  
  tad$gene.count <- 0
  tad$genes <- NA
  for(i in 1:nrow(tad)) {
    g <- genes.data[genes.data$chr == tad$chr[i], ] 
    hits <- g[g$tss > tad$bp.start[i] & g$tss < tad$bp.stop[i], ]
    tad$gene.count[i] <- nrow(hits)
    tad$genes[i] <- paste(hits$Gene.name, collapse = ';')
  }
  
  tad$range <- paste(tad$bp.start, tad$bp.stop, sep = '-')
  tad$ID <- paste(tad$chr, 1:nrow(tad), tad$range, sep = '.')
  
  tad$length <- tad$bp.stop-tad$bp.start
  tad$genes <- paste(tad$genes, ";", sep = '')
  return (tad)
}


# Description:
# Gene is in TAD if it overlaps (even partially) the boundary of a TAD
getTadStats <- function(tad.dir, genes.data) {
  # remove duplicated genes
  genes.data <- genes.data[duplicated(genes.data$Gene.name) == F, ]
  
  tad <- read.csv(file = tad.dir, sep = '\t', header = F, as.is = T)
  colnames(tad) <- c("chr", "bp.start", "bp.stop")
  tad$mb.length <- (tad$bp.stop-tad$bp.start)/10^6
  tad$chr <- factor(x = tad$chr, levels = getChrs())
  tad$range <- paste(tad$bp.start, tad$bp.stop, sep = '-')
  tad$ID <- paste(tad$chr, 1:nrow(tad), tad$range, sep = '.')
  tad$gene.count <- 0
  tad$de.gene.count <- 0
  tad$genes <- 0
  
  
  tad.stats <- c()
  for(chr in unique(tad$chr)) {
    temp.tad <- tad[tad$chr == chr, ]
    temp.genes <- genes.data[genes.data$chr == chr, ]
    
    q <- IRanges(start = temp.tad$bp.start, end = temp.tad$bp.stop)
    s <- IRanges(start = temp.genes$bp.start, end = temp.genes$bp.stop)
    m <- findOverlaps(query = q, subject = s, type = "any")
    m <- data.frame(as.matrix(m))
    if(nrow(m) != 0) {
      m$ID <- temp.tad$ID[m$queryHits]
      m$genes <- temp.genes$Gene.name[m$subjectHits]
      if("significant" %in% colnames(temp.genes) == TRUE) {
        m$significant <- temp.genes$significant[m$subjectHits]
      }
      m$gene.count <- 1
      m.tad.counts <- aggregate(gene.count~ID, data = m, FUN = sum)
      m.tad.de.counts <- aggregate(significant~ID, data = m, FUN = sum)
      m.tad.genes <- aggregate(genes~ID, data = m, FUN = base::paste, collapse = ';')
      m <- merge(x = m.tad.genes, y = m.tad.counts, by = "ID")
      m <- merge(x = m, y = m.tad.de.counts, by = "ID")
      tad.stats <- rbind(tad.stats, m)
    }
    else {
      m <- data.frame(ID = temp.tad$ID, 
                      genes = NA, 
                      gene.count = NA, 
                      significant = NA)
      tad.stats <- rbind(tad.stats, m)
    }
  }
  
  key <- do.call(rbind, strsplit(x = tad.stats$ID, split = '\\.'))
  key.range <- do.call(rbind, strsplit(x = key[, 3], split = '\\-'))
  tad.stats$chr <- key[, 1]
  tad.stats$bp.start <- as.numeric(key.range[, 1])
  tad.stats$bp.stop <- as.numeric(key.range[, 2])
  return (tad.stats)
}






# Description:
# Get formatted TAD data
getTadData <- function(tad.dir) {
  tad <- read.csv(file = tad.dir, sep = '\t', header = F, as.is = T)
  colnames(tad) <- c("chr", "bp.start", "bp.stop")
  tad$bp.mid <- tad$bp.start+(tad$bp.stop-tad$bp.start)/2
  tad$mb.length <- (tad$bp.stop-tad$bp.start)/10^6
  tad$chr <- factor(x = tad$chr, levels = getChrs())
  tad$range <- paste(tad$bp.start, tad$bp.stop, sep = '-')
  tad$ID <- paste(tad$chr, 1:nrow(tad), tad$range, sep = '.')
  
  return (tad)
}

