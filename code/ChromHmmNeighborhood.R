#############################
## Neighborhood enrichment ##
#############################

# 0. get utils
source("code/util/ChromHmmUtilities.R")
source("code/util/GeneralUtilities.R")


# Description:
# - bin/state-centrics analysis
# - computes different kind of overlaps (state-specific/region-specific/ChromHMM enrichment-fold)
# - regional (external annotation) data taken from ChromHMM
getOverlapData <- function(states) {
  # c. create overlap
  states$key <- paste(states$chr, states$bp.start, sep = '.')
  
  bed.files <- list.files(path = "output/chromhmm_neighborhood/", 
                          pattern = "mm10.bed|fantom_", full.names = T)
  short.bed.files <- list.files(path = "output/chromhmm_neighborhood/", 
                                pattern = "mm10.bed|fantom_", full.names = F)
  for(i in 1:length(bed.files)) {
    file <- read.csv(file = bed.files[i], sep = "\t", header = F, as.is = T)
    file$key <- paste(file$V1, file$V2, sep = '.')
    states[, short.bed.files[i]] <- 0
    states[states$key %in% file$key, short.bed.files[i]] <- 1
    cat(i, "\n")
  }
  return (states)
}



# 1. combine bins and hiddenDomain peaks
header <- get(load(file = "output/header.RData"))



# a. create header
write.table(x = header, file = "output/chromhmm_neighborhood/header.bed", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = F)




# b. ZNF genes
genes <- getGenes(dir = "input/genes_biomart/genes.full.txt")
znf.genes <- genes[which(regexpr(pattern = "Zfp", text = genes$Gene.name, ignore.case = T) != -1), ]
znf.genes <- genes[which(regexpr(pattern = "zinc finger protein", 
                                 text = genes$Gene.description, 
                                 ignore.case = T) != -1), ]

write.table(x = znf.genes[, c("chr", "Gene.start..bp.", "Gene.end..bp.")], 
            file = "input/mm10_annotation/Zfp.mm10.bed", sep = "\t", 
            row.names = F, col.names = F, quote = F)




# c. create bed file intersection 
files <- list.files(path = "input/mm10_annotation/", pattern = "bed", full.names = T)
short.files <- list.files(path = "input/mm10_annotation/", pattern = "bed", full.names = F)
for(i in 1:length(files)) {
  cat(i, "\n")
  system(command = paste("bedtools intersect -a output/chromhmm_neighborhood/header.bed -b ", 
                         files[i], " -u > output/chromhmm_neighborhood/", short.files[i], sep = ''))
}



# optional: d. fantom file intersection
files <- list.files(path = "input/fantom/", pattern = "bed", full.names = T)
short.files <- list.files(path = "input/fantom/", pattern = "bed", full.names = F)
for(i in 1:length(files)) {
  cat(i, "\n")
  system(command = paste("bedtools intersect -a output/chromhmm_neighborhood/header.bed -b ", 
                         files[i], " -u > output/chromhmm_neighborhood/", short.files[i], sep = ''))
}
rm(i, files, short.files)




# 2. get overlap data
state.data <- get(load(file = "output/states.15.RData"))
overlap <- getOverlapData(states = state.data)
save(overlap, file = "output/chromhmm_neighborhood/overlap.15.RData")
rm(state.data, overlap)






# A = # of bases in state
# B = # of bases in external annotation
# C = # of bases in state AND external annotation
# D = # of bases in genome
# p = fold-enrichment = (C/A)/(B/D)
# 3. replicate ChromHMM enrichment
overlap <- get(load(file = "output/chromhmm_neighborhood/overlap.15.RData"))
regions <- colnames(overlap)[regexpr(
  pattern = "mm10|fantom_forelimb|fantom_allEnhancers", 
  text = colnames(overlap)) != -1]


states <- 1:15
chrom.stats <- c()
for(j in 1:length(regions)) {
  B <- sum(overlap[, regions[j]])
  for(i in 1:length(states)) {
    p.sum <- data.frame(A = sum(overlap$states.P == states[i]), B = B, D = nrow(overlap),
                        C = sum(overlap$states.P == states[i] & overlap[, regions[j]] == 1),
                        genomic.region = regions[j], cell.type = "PC", state = states[i], 
                        stringsAsFactors = F)
    h.sum <- data.frame(A = sum(overlap$states.H == states[i]), B = B, D = nrow(overlap),
                        C = sum(overlap$states.H == states[i] & overlap[, regions[j]] == 1),
                        genomic.region = regions[j], cell.type = "HC", state = states[i], 
                        stringsAsFactors = F)
    chrom.stats <- rbind(chrom.stats, p.sum)
    chrom.stats <- rbind(chrom.stats, h.sum)
  }
  cat(j, "\n")
}

chrom.stats$p <- (chrom.stats$C/chrom.stats$A)/(chrom.stats$B/chrom.stats$D)
chrom.stats$genomic.region <- gsub(pattern = ".mm10.bed.gz|.mm10.bed|_E14.bed|.bed", 
                                   replacement = '', x = chrom.stats$genomic.region)
chrom.stats$genomic.region[chrom.stats$genomic.region == "enhancers"] <- "Vista"
chrom.stats$genomic.region[chrom.stats$genomic.region == "fantom_allEnhancers"] <- "Fantom all"
chrom.stats$genomic.region[chrom.stats$genomic.region == "fantom_forelimb_embryo"] <- "Fantom forelimb"

chrom.stats$genomic.region <- factor(x = chrom.stats$genomic.region,
                                     levels = c("RefSeqGene", "RefSeqExon", 
                                                "RefSeqTSS", "RefSeqTSS2kb", 
                                                "RefSeqTES", "CpGIsland",
                                                "Vista", "Fantom all", 
                                                "Fantom forelimb", "Zfp"))
      
chrom.stats$cell.type <- factor(x = chrom.stats$cell.type, levels = c("PC", "HC"))

g <- ggplot(data = chrom.stats)+
  facet_wrap(facets = ~cell.type, ncol = 2)+
  geom_tile(aes(x = genomic.region, y = as.factor(state), fill = p), col = "black")+
  geom_text(aes(x = genomic.region, y = as.factor(state), label = round(x = p, digits = 2)), 
            col = "black", size = 3.55)+
  theme_bw(base_size = 12)+
  scale_fill_continuous(name = "Fold enrichment", low = "white", high = "#4682b4")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        legend.position = "top")+
  ylab(label = "State")+
  xlab(label = "")

ggsave(filename = "output/chromhmm_plots/Enrichment.15.tiff", plot = g,
       device = "tiff", width = 8, height = 7, dpi = 300)








# 3b
state.data <- get(load(file = "output/states.LB.15.RData"))
overlap <- getOverlapData(states = state.data)
save(overlap, file = "output/chromhmm_neighborhood/overlap.LB.15.RData")
rm(state.data, overlap)


# 4b. replicate ChromHMM enrichment (LB)
overlap <- get(load(file = "output/chromhmm_neighborhood/overlap.LB.15.RData"))
regions <- colnames(overlap)[regexpr(
  pattern = "mm10|fantom_forelimb|fantom_allEnhancers", 
  text = colnames(overlap)) != -1]

states <- 1:15
chrom.stats <- c()
for(j in 1:length(regions)) {
  B <- sum(overlap[, regions[j]])
  for(i in 1:length(states)) {
    p.sum <- data.frame(A = sum(overlap$states.P == states[i]), B = B, D = nrow(overlap),
                        C = sum(overlap$states.P == states[i] & overlap[, regions[j]] == 1),
                        genomic.region = regions[j], cell.type = "PC", state = states[i], stringsAsFactors = F)
    h.sum <- data.frame(A = sum(overlap$states.H == states[i]), B = B, D = nrow(overlap),
                        C = sum(overlap$states.H == states[i] & overlap[, regions[j]] == 1),
                        genomic.region = regions[j], cell.type = "HC", state = states[i], stringsAsFactors = F)
    chrom.stats <- rbind(chrom.stats, p.sum)
    chrom.stats <- rbind(chrom.stats, h.sum)
  }
  cat(j, "\n")
}

chrom.stats$p <- (chrom.stats$C/chrom.stats$A)/(chrom.stats$B/chrom.stats$D)
chrom.stats$genomic.region <- gsub(pattern = ".mm10.bed.gz|.mm10.bed|_E14.bed|.bed", 
                                   replacement = '', x = chrom.stats$genomic.region)
chrom.stats$genomic.region[chrom.stats$genomic.region == "enhancers"] <- "Vista"
chrom.stats$genomic.region[chrom.stats$genomic.region == "fantom_allEnhancers"] <- "Fantom all"
chrom.stats$genomic.region[chrom.stats$genomic.region == "fantom_forelimb_embryo"] <- "Fantom forelimb"

chrom.stats$genomic.region <- factor(x = chrom.stats$genomic.region,
                                     levels = c("RefSeqGene", "RefSeqExon", 
                                                "RefSeqTSS", "RefSeqTSS2kb", 
                                                "RefSeqTES", "CpGIsland",
                                                "Vista", "Fantom all", 
                                                "Fantom forelimb", "Zfp"))

chrom.stats$cell.type <- factor(x = chrom.stats$cell.type, levels = c("PC", "HC"))

g <- ggplot(data = chrom.stats)+
  facet_wrap(facets = ~cell.type, ncol = 2)+
  geom_tile(aes(x = genomic.region, y = as.factor(state), fill = p), col = "black")+
  geom_text(aes(x = genomic.region, y = as.factor(state), label = round(x = p, digits = 2)), 
            col = "black", size = 3.55)+
  theme_bw(base_size = 12)+
  scale_fill_continuous(name = "Fold enrichment", low = "white", high = "#4682b4")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        legend.position = "top")+
  ylab(label = "State")+
  xlab(label = "")


ggsave(filename = "output/chromhmm_plots/Enrichment.LB.15.tiff", plot = g,
       device = "tiff", width = 8, height = 7, dpi = 300)

