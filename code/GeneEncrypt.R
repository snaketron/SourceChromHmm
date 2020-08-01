##################################
## Map ChromHMM states to genes ##
##################################


# require code/libs
source("code/util/GeneralUtilities.R")
require("IRanges")
require("parallel")



# Description:
# Encode island presence information for each gene in 3 different genomic regions
getGcsc3 <- function(chr,
                     island.data, 
                     gene.data, 
                     tss.buffer = 2500, 
                     enhancer.buffer = 2*10^4) {
  
  
  # subset by chr
  island.data <- island.data[island.data$chr == chr, ]
  gene.data <- gene.data[gene.data$chr == chr, ]
  
  
  # ehancers
  gene.data$bp.start.enhancer <- gene.data$bp.start-enhancer.buffer
  gene.data$bp.stop.enhancer <- gene.data$bp.stop+enhancer.buffer
  q <- IRanges(start = gene.data$bp.start.enhancer, end = gene.data$bp.stop.enhancer)
  s <- IRanges(start = island.data$bp.start, end = island.data$bp.stop)
  m <- findOverlaps(query = q, subject = s, type = "any")
  m <- data.frame(as.matrix(m))
  m$state <- island.data$state[m$subjectHits]
  m$gene <- gene.data$Gene.name[m$queryHits]
  m$gene.id <- gene.data$Gene.stable.ID[m$queryHits]
  m$strand <- gene.data$Strand[m$queryHits]
  m$chr <- gene.data$chr[m$queryHits]
  m$count <- 1
  m.enhancer <- aggregate(count~chr+gene.id+strand+state, data = m, FUN = sum)
  m.enhancer$region <- "enhancer"
  
  
  # tss
  gene.data$bp.start.tss <- gene.data$bp.start-tss.buffer
  gene.data$bp.stop.tss <- gene.data$bp.stop+tss.buffer
  q <- IRanges(start = gene.data$bp.start.tss, end = gene.data$bp.stop.tss)
  s <- IRanges(start = island.data$bp.start, end = island.data$bp.stop)
  m <- findOverlaps(query = q, subject = s, type = "any")
  m <- data.frame(as.matrix(m))
  m$state <- island.data$state[m$subjectHits]
  m$gene <- gene.data$Gene.name[m$queryHits]
  m$gene.id <- gene.data$Gene.stable.ID[m$queryHits]
  m$strand <- gene.data$Strand[m$queryHits]
  m$chr <- gene.data$chr[m$queryHits]
  m$count <- 1
  m.tss <- aggregate(count~chr+gene.id+strand+state, data = m, FUN = sum)
  m.tss$region <- "tss"
  
  
  # gene
  q <- IRanges(start = gene.data$bp.start, end = gene.data$bp.stop)
  s <- IRanges(start = island.data$bp.start, end = island.data$bp.stop)
  m <- findOverlaps(query = q, subject = s, type = "any")
  m <- data.frame(as.matrix(m))
  m$state <- island.data$state[m$subjectHits]
  m$gene <- gene.data$Gene.name[m$queryHits]
  m$gene.id <- gene.data$Gene.stable.ID[m$queryHits]
  m$strand <- gene.data$Strand[m$queryHits]
  m$chr <- gene.data$chr[m$queryHits]
  m$count <- 1
  m.gene <- aggregate(count~chr+gene.id+strand+state, data = m, FUN = sum)
  m.gene$region <- "gene"
  
  
  return (rbind(m.enhancer, m.tss, m.gene))
}


# gene data
gene.data <- getGenes(dir = "input/genes_biomart/genes.full.txt")



#### mutations ####
island.data <- get(load(file = "output/chromhmm_islands/mutation.islands.15.RData"))
rm(islands)

m <- mclapply(X = unique(island.data$chr), 
              FUN = getGcsc3, 
              island.data = island.data, 
              gene.data = gene.data, 
              tss.buffer = 2500, 
              enhancer.buffer = 2*10^4, 
              mc.cores = 20)
m <- do.call(rbind, m)
save(m, file = "output/gene_centric/gcsc3.mutations.RData")
cat("Done \n")
rm(gene.data, island.data, chrs)






#### states ####


# gene data
gene.data <- getGenes(dir = "input/genes_biomart/genes.full.txt")


# PC
pc.islands <- get(load(file = "output/chromhmm_islands/PC.islands.15.RData"))
rm(islands.PC)

pc <- mclapply(X = unique(pc.islands$chr), 
              FUN = getGcsc3, 
              island.data = pc.islands, 
              gene.data = gene.data, 
              tss.buffer = 2500, 
              enhancer.buffer = 2*10^4, 
              mc.cores = 20)
pc <- do.call(rbind, pc)
pc$cell.type = "PC"


# HC
hc.islands <- get(load(file = "output/chromhmm_islands/HC.islands.15.RData"))
rm(islands.HC)

hc <- mclapply(X = unique(hc.islands$chr), 
               FUN = getGcsc3, 
               island.data = hc.islands, 
               gene.data = gene.data, 
               tss.buffer = 2500, 
               enhancer.buffer = 2*10^4, 
               mc.cores = 20)
hc <- do.call(rbind, hc)
hc$cell.type = "HC"


o <- rbind(pc, hc)
save(o, file = "output/gene_centric/gcsc3.states.RData")
cat("Done \n")
rm(pc, hc, pc.islands, hc.islands)
