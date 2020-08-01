# TadDelta computation

source("code/util/GeneralUtilities.R")
source("code/util/ChromHmmUtilities.R")
source("code/util/TadStats.R")
source("code/util/EnhancerUtillities.R")
require(ggplot2)
require(ggrepel)


# H3K27ac peaks not overlapping TSSs
ePC <- get(load(file = "output/enhancers/ePC.RData"))
eHC <- get(load(file = "output/enhancers/eHC.RData"))
ePC <- ePC[ePC$overlap.TSS == 0, ]
eHC <- eHC[eHC$overlap.TSS == 0, ]


# genes
gene.data <- getGenes(dir = "input/genes_biomart/genes.full.txt")
gene.data <- gene.data[gene.data$Gene.type == "protein_coding", ]


# TADs
tad.hind3 <- read.csv(file = "input/tad/total.HindIII.combined.domain",
                      header = FALSE, sep = "\t", as.is = TRUE)
colnames(tad.hind3) <- c("chr", "bp.start", "bp.stop")
tad.hind3$TAD.ID <- "total.HindIII.combined.domain"
tad.hind3$id <- "TAD: HindIII"

tad.ncoi <- read.csv(file = "input/tad/total.NcoI.domains",
                     header = FALSE, sep = "\t", as.is = TRUE)
colnames(tad.ncoi) <- c("chr", "bp.start", "bp.stop")
tad.ncoi$TAD.ID <- "total.NcoI.domains"
tad.ncoi$id <- "TAD: NcoI"

tads <- rbind(tad.hind3, tad.ncoi)
rm(tad.hind3, tad.ncoi)

tads$count.PC <- 0
tads$count.HC <- 0
tads$genes <- ''
for(i in 1:nrow(tads)) {
  tads$count.PC[i] <- sum((ePC$chr == tads$chr[i] & 
                             ePC$bp.mid >= tads$bp.start[i] & 
                             ePC$bp.mid <= tads$bp.stop[i]) |
                            (ePC$chr == tads$chr[i] & 
                               ePC$bp.start >= tads$bp.start[i] & 
                               ePC$bp.start <= tads$bp.stop[i]) |
                            (ePC$chr == tads$chr[i] & 
                               ePC$bp.stop >= tads$bp.start[i] & 
                               ePC$bp.stop <= tads$bp.stop[i]))
  tads$count.HC[i] <- sum((eHC$chr == tads$chr[i] & 
                             eHC$bp.mid >= tads$bp.start[i] & 
                             eHC$bp.mid <= tads$bp.stop[i]) |
                            (eHC$chr == tads$chr[i] & 
                               eHC$bp.start >= tads$bp.start[i] & 
                               eHC$bp.start <= tads$bp.stop[i]) |
                            (eHC$chr == tads$chr[i] & 
                               eHC$bp.stop >= tads$bp.start[i] & 
                               eHC$bp.stop <= tads$bp.stop[i]))
  
  g <- which(gene.data$chr == tads$chr[i] & 
               gene.data$tss >= tads$bp.start[i] & 
               gene.data$tss <= tads$bp.stop[i])
  if(length(g) != 0) {
    tads$genes[i] <- paste(gene.data$Gene.name[g], collapse = ',')
  }
}
tads$delta <- tads$count.PC-tads$count.HC
tads$lfc <- log(tads$count.PC/tads$count.HC)
rm(g, i)




# create ranks
x <- tads[tads$id == "TAD: HindIII", ]
x <- x[order(x$delta, decreasing = T), ]
x$rank.PC <- 1:nrow(x)
x <- x[order(x$delta, decreasing = F), ]
x$rank.HC <- 1:nrow(x)


y <- tads[tads$id == "TAD: NcoI", ]
y <- y[order(y$delta, decreasing = TRUE), ]
y$rank.PC <- 1:nrow(y)
y <- y[order(y$delta, decreasing = F), ]
y$rank.HC <- 1:nrow(y)


tads <- rbind(x, y)
rm(x, y)
tads$bit <- ifelse(test = tads$rank.PC <= 5 | tads$rank.HC <= 5, 
                   yes = T, no = F)


tads$genes.wrap <- gsub(pattern = ',', replacement = '\n', x = tads$genes)




# plots
g1 <- ggplot(data = tads)+
  facet_wrap(facets = ~id, nrow = 1)+
  geom_point(aes(x = count.PC, y = count.HC,  col = bit))+
  geom_density_2d(aes(x = count.PC, y = count.HC))+
  geom_abline(slope = 1, intercept = 0, col = "red")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_manual(values = c("black", "red"))+
  xlab(label = "# of PC H3K27ac peaks (non-TSS) in TAD")+
  ylab(label = "# of HC H3K27ac peaks (non-TSS) in TAD")


g2 <- ggplot(data = tads)+
  facet_wrap(facets = ~id, nrow = 1, labeller = label_parsed)+
  geom_point(aes(x = rank.PC, y = delta, col = bit), size = 1)+
  geom_label_repel(data = tads[tads$bit==T,], size = 2, fill = NA, 
                   min.segment.length = 0, force = 2, max.iter = 5000,
                  aes(x = rank.PC, y = delta, label = genes.wrap))+
  theme_bw()+
  theme(legend.position = "none")+
  xlab(label = "Ranked TADs")+
  scale_color_manual(values = c("darkgray", "red"))+
  ylab(label = "Delta TAD peaks (PC-HC)")


g <- gridExtra::arrangeGrob(g1, g2, heights = c(2, 3))
plot(g)


ggsave(filename = paste0("output/enhancers/tad_delta.tiff"),
       plot = g, device = "tiff", width = 10, height = 10, 
       dpi = 300)

save(tads, file = "output/enhancers/tad_delta.RData")


# tables
tads$bit <- NULL
tads$lfc <- NULL
tads$rank.HC <- NULL
tads$rank.PC <- NULL
tads$genes.wrap <- NULL
tads$TAD.ID <- NULL
write.table(x = tads, file = "output/enhancers/tad_delta.tsv",
            sep = '\t', row.names = F, quote = F)
