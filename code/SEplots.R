# Plots for figures involving SEs and TADs

source("code/util/GeneralUtilities.R")
source("code/util/ChromHmmUtilities.R")
source("code/util/TadStats.R")
source("code/util/EnhancerUtillities.R")
require("ggplot2")
require("gridExtra")
require("reshape2")
require("tidyr")
require("ggplot2")
require("ggrepel")
require("gridExtra")


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



# genes: Sox5, Sox9
gene.data <- getGenes(dir = "input/genes_biomart/genes.full.txt")
gene.data <- gene.data[gene.data$Gene.type == "protein_coding", ]
sox5 <- gene.data[gene.data$Gene.name == "Sox5", ]
sox9 <- gene.data[gene.data$Gene.name == "Sox9", ]
col10 <- gene.data[gene.data$Gene.name == "Col10a1", ]



# get overlapping TADs
tads.sox5 <- tads[tads$bp.start <= sox5$tss & 
                    tads$bp.stop >= sox5$tss & 
                    tads$chr == sox5$chr, ] 

tads.sox9 <- tads[tads$bp.start <= sox9$tss & 
                    tads$bp.stop >= sox9$tss & 
                    tads$chr == sox9$chr, ] 

tads.col10 <- tads[tads$bp.start <= col10$tss & 
                     tads$bp.stop >= col10$tss & 
                     tads$chr == col10$chr, ] 
rm(tads)


# get genes with TSS within either TAD
genes.sox5 <- gene.data[which((gene.data$tss >= tads.sox5$bp.start[1] & 
                                 gene.data$tss <= tads.sox5$bp.stop[1] & 
                                 gene.data$chr == tads.sox5$chr[1]) |
                                (gene.data$bp.start >= tads.sox5$bp.start[1] & 
                                   gene.data$bp.start <= tads.sox5$bp.stop[1] & 
                                   gene.data$chr == tads.sox5$chr[1]) |
                                (gene.data$bp.stop >= tads.sox5$bp.start[1] & 
                                   gene.data$bp.stop <= tads.sox5$bp.stop[1] & 
                                   gene.data$chr == tads.sox5$chr[1]) |
                                (gene.data$tss >= tads.sox5$bp.start[2] & 
                                   gene.data$tss <= tads.sox5$bp.stop[2] & 
                                   gene.data$chr == tads.sox5$chr[2]) |
                                (gene.data$bp.start >= tads.sox5$bp.start[2] & 
                                   gene.data$bp.start <= tads.sox5$bp.stop[2] & 
                                   gene.data$chr == tads.sox5$chr[2]) |
                                (gene.data$bp.stop >= tads.sox5$bp.start[2] & 
                                   gene.data$bp.stop <= tads.sox5$bp.stop[2] & 
                                   gene.data$chr == tads.sox5$chr[2])), ] 


genes.sox9 <- gene.data[which((gene.data$tss >= tads.sox9$bp.start[1] & 
                                 gene.data$tss <= tads.sox9$bp.stop[1] & 
                                 gene.data$chr == tads.sox9$chr[1]) |
                                (gene.data$bp.start >= tads.sox9$bp.start[1] & 
                                   gene.data$bp.start <= tads.sox9$bp.stop[1] & 
                                   gene.data$chr == tads.sox9$chr[1]) |
                                (gene.data$bp.stop >= tads.sox9$bp.start[1] & 
                                   gene.data$bp.stop <= tads.sox9$bp.stop[1] & 
                                   gene.data$chr == tads.sox9$chr[1]) |
                                (gene.data$tss >= tads.sox9$bp.start[2] & 
                                   gene.data$tss <= tads.sox9$bp.stop[2] & 
                                   gene.data$chr == tads.sox9$chr[2]) |
                                (gene.data$bp.start >= tads.sox9$bp.start[2] & 
                                   gene.data$bp.start <= tads.sox9$bp.stop[2] & 
                                   gene.data$chr == tads.sox9$chr[2]) |
                                (gene.data$bp.stop >= tads.sox9$bp.start[2] & 
                                   gene.data$bp.stop <= tads.sox9$bp.stop[2] & 
                                   gene.data$chr == tads.sox9$chr[2])), ] 


genes.col10 <- gene.data[which((gene.data$tss >= tads.col10$bp.start[1] & 
                                  gene.data$tss <= tads.col10$bp.stop[1] & 
                                  gene.data$chr == tads.col10$chr[1]) |
                                 (gene.data$bp.start >= tads.col10$bp.start[1] & 
                                    gene.data$bp.start <= tads.col10$bp.stop[1] & 
                                    gene.data$chr == tads.col10$chr[1]) |
                                 (gene.data$bp.stop >= tads.col10$bp.start[1] & 
                                    gene.data$bp.stop <= tads.col10$bp.stop[1] & 
                                    gene.data$chr == tads.col10$chr[1]) |
                                 (gene.data$tss >= tads.col10$bp.start[2] & 
                                    gene.data$tss <= tads.col10$bp.stop[2] & 
                                    gene.data$chr == tads.col10$chr[2]) |
                                 (gene.data$bp.start >= tads.col10$bp.start[2] & 
                                    gene.data$bp.start <= tads.col10$bp.stop[2] & 
                                    gene.data$chr == tads.col10$chr[2]) |
                                 (gene.data$bp.stop >= tads.col10$bp.start[2] & 
                                    gene.data$bp.stop <= tads.col10$bp.stop[2] & 
                                    gene.data$chr == tads.col10$chr[2])), ] 
rm(gene.data)



# Vista: get vista elements with midpoints within either TAD
vista <- get(load(file = "input/vista/vista.RData"))
vista$bp.mid <- vista$bp.start + (vista$bp.stop-vista$bp.start)/2

vista.sox5 <- vista[which((vista$bp.mid >= tads.sox5$bp.start[1] & 
                             vista$bp.mid <= tads.sox5$bp.stop[1] & 
                             vista$chr == tads.sox5$chr[1]) |
                            (vista$bp.start >= tads.sox5$bp.start[1] & 
                               vista$bp.start <= tads.sox5$bp.stop[1] & 
                               vista$chr == tads.sox5$chr[1]) |
                            (vista$bp.stop >= tads.sox5$bp.start[1] & 
                               vista$bp.stop <= tads.sox5$bp.stop[1] & 
                               vista$chr == tads.sox5$chr[1]) |
                            (vista$bp.mid >= tads.sox5$bp.start[2] & 
                               vista$bp.mid <= tads.sox5$bp.stop[2] & 
                               vista$chr == tads.sox5$chr[2]) |
                            (vista$bp.start >= tads.sox5$bp.start[2] & 
                               vista$bp.start <= tads.sox5$bp.stop[2] & 
                               vista$chr == tads.sox5$chr[2]) |
                            (vista$bp.stop >= tads.sox5$bp.start[2] & 
                               vista$bp.stop <= tads.sox5$bp.stop[2] & 
                               vista$chr == tads.sox5$chr[2])), ] 

vista.sox9 <- vista[which((vista$bp.mid >= tads.sox9$bp.start[1] & 
                             vista$bp.mid <= tads.sox9$bp.stop[1] & 
                             vista$chr == tads.sox9$chr[1]) |
                            (vista$bp.start >= tads.sox9$bp.start[1] & 
                               vista$bp.start <= tads.sox9$bp.stop[1] & 
                               vista$chr == tads.sox9$chr[1]) |
                            (vista$bp.stop >= tads.sox9$bp.start[1] & 
                               vista$bp.stop <= tads.sox9$bp.stop[1] & 
                               vista$chr == tads.sox9$chr[1]) |
                            (vista$bp.mid >= tads.sox9$bp.start[2] & 
                               vista$bp.mid <= tads.sox9$bp.stop[2] & 
                               vista$chr == tads.sox9$chr[2]) |
                            (vista$bp.start >= tads.sox9$bp.start[2] & 
                               vista$bp.start <= tads.sox9$bp.stop[2] & 
                               vista$chr == tads.sox9$chr[2]) |
                            (vista$bp.stop >= tads.sox9$bp.start[2] & 
                               vista$bp.stop <= tads.sox9$bp.stop[2] & 
                               vista$chr == tads.sox9$chr[2])), ] 

vista.col10 <- vista[which((vista$bp.mid >= tads.col10$bp.start[1] & 
                              vista$bp.mid <= tads.col10$bp.stop[1] & 
                              vista$chr == tads.col10$chr[1]) |
                             (vista$bp.start >= tads.col10$bp.start[1] & 
                                vista$bp.start <= tads.col10$bp.stop[1] & 
                                vista$chr == tads.col10$chr[1]) |
                             (vista$bp.stop >= tads.col10$bp.start[1] & 
                                vista$bp.stop <= tads.col10$bp.stop[1] & 
                                vista$chr == tads.col10$chr[1]) |
                             (vista$bp.mid >= tads.col10$bp.start[2] & 
                                vista$bp.mid <= tads.col10$bp.stop[2] & 
                                vista$chr == tads.col10$chr[2]) |
                             (vista$bp.start >= tads.col10$bp.start[2] & 
                                vista$bp.start <= tads.col10$bp.stop[2] & 
                                vista$chr == tads.col10$chr[2]) |
                             (vista$bp.stop >= tads.col10$bp.start[2] & 
                                vista$bp.stop <= tads.col10$bp.stop[2] & 
                                vista$chr == tads.col10$chr[2])), ] 
rm(vista)





# SEs: get SE elements with midpoints within either TAD
SE.data <- get(load(file = "output/enhancers/SE.RData"))
SE.data$bp.mid <- SE.data$bp.start + (SE.data$bp.stop-SE.data$bp.start)/2
SE.data <- SE.data[SE.data$count >= 5, ]

SE.sox5 <- SE.data[which((SE.data$bp.mid >= tads.sox5$bp.start[1] & 
                            SE.data$bp.mid <= tads.sox5$bp.stop[1] & 
                            SE.data$chr == tads.sox5$chr[1]) |
                           (SE.data$bp.start >= tads.sox5$bp.start[1] & 
                              SE.data$bp.start <= tads.sox5$bp.stop[1] & 
                              SE.data$chr == tads.sox5$chr[1]) |
                           (SE.data$bp.stop >= tads.sox5$bp.start[1] & 
                              SE.data$bp.stop <= tads.sox5$bp.stop[1] & 
                              SE.data$chr == tads.sox5$chr[1]) |
                           (SE.data$bp.start >= tads.sox5$bp.start[2] & 
                              SE.data$bp.start <= tads.sox5$bp.stop[2] & 
                              SE.data$chr == tads.sox5$chr[2]) |
                           (SE.data$bp.stop >= tads.sox5$bp.start[2] & 
                              SE.data$bp.stop <= tads.sox5$bp.stop[2] & 
                              SE.data$chr == tads.sox5$chr[2]) |
                           (SE.data$bp.mid >= tads.sox5$bp.start[2] & 
                              SE.data$bp.mid <= tads.sox5$bp.stop[2] & 
                              SE.data$chr == tads.sox5$chr[2])), ] 

SE.sox9 <- SE.data[which((SE.data$bp.mid >= tads.sox9$bp.start[1] & 
                            SE.data$bp.mid <= tads.sox9$bp.stop[1] & 
                            SE.data$chr == tads.sox9$chr[1]) |
                           (SE.data$bp.start >= tads.sox9$bp.start[1] & 
                              SE.data$bp.start <= tads.sox9$bp.stop[1] & 
                              SE.data$chr == tads.sox9$chr[1]) |
                           (SE.data$bp.stop >= tads.sox9$bp.start[1] & 
                              SE.data$bp.stop <= tads.sox9$bp.stop[1] & 
                              SE.data$chr == tads.sox9$chr[1]) |
                           (SE.data$bp.start >= tads.sox9$bp.start[2] & 
                              SE.data$bp.start <= tads.sox9$bp.stop[2] & 
                              SE.data$chr == tads.sox9$chr[2]) |
                           (SE.data$bp.stop >= tads.sox9$bp.start[2] & 
                              SE.data$bp.stop <= tads.sox9$bp.stop[2] & 
                              SE.data$chr == tads.sox9$chr[2]) |
                           (SE.data$bp.mid >= tads.sox9$bp.start[2] & 
                              SE.data$bp.mid <= tads.sox9$bp.stop[2] & 
                              SE.data$chr == tads.sox9$chr[2])), ]

SE.col10 <- SE.data[which((SE.data$bp.mid >= tads.col10$bp.start[1] & 
                            SE.data$bp.mid <= tads.col10$bp.stop[1] & 
                            SE.data$chr == tads.col10$chr[1]) |
                           (SE.data$bp.start >= tads.col10$bp.start[1] & 
                              SE.data$bp.start <= tads.col10$bp.stop[1] & 
                              SE.data$chr == tads.col10$chr[1]) |
                           (SE.data$bp.stop >= tads.col10$bp.start[1] & 
                              SE.data$bp.stop <= tads.col10$bp.stop[1] & 
                              SE.data$chr == tads.col10$chr[1]) |
                           (SE.data$bp.start >= tads.col10$bp.start[2] & 
                              SE.data$bp.start <= tads.col10$bp.stop[2] & 
                              SE.data$chr == tads.col10$chr[2]) |
                           (SE.data$bp.stop >= tads.col10$bp.start[2] & 
                              SE.data$bp.stop <= tads.col10$bp.stop[2] & 
                              SE.data$chr == tads.col10$chr[2]) |
                           (SE.data$bp.mid >= tads.col10$bp.start[2] & 
                              SE.data$bp.mid <= tads.col10$bp.stop[2] & 
                              SE.data$chr == tads.col10$chr[2])), ]
rm(SE.data)



# Peaks: get H3K27ac peaks with midpoints within either TAD
ePC <- get(load(file = "output/enhancers/ePC.RData"))
eHC <- get(load(file = "output/enhancers/eHC.RData"))
e <- rbind(ePC, eHC)
e <- e[e$overlap.TSS == 0, ]
rm(ePC, eHC)


e.sox5 <- e[which((e$bp.mid >= tads.sox5$bp.start[1] & 
                     e$bp.mid <= tads.sox5$bp.stop[1] & 
                     e$chr == tads.sox5$chr[1]) |
                    (e$bp.start >= tads.sox5$bp.start[1] & 
                       e$bp.start <= tads.sox5$bp.stop[1] & 
                       e$chr == tads.sox5$chr[1]) |
                    (e$bp.stop >= tads.sox5$bp.start[1] & 
                       e$bp.stop <= tads.sox5$bp.stop[1] & 
                       e$chr == tads.sox5$chr[1]) |
                    (e$bp.start >= tads.sox5$bp.start[2] & 
                       e$bp.start <= tads.sox5$bp.stop[2] & 
                       e$chr == tads.sox5$chr[2]) |
                    (e$bp.stop >= tads.sox5$bp.start[2] & 
                       e$bp.stop <= tads.sox5$bp.stop[2] & 
                       e$chr == tads.sox5$chr[2]) |
                    (e$bp.mid >= tads.sox5$bp.start[2] & 
                       e$bp.mid <= tads.sox5$bp.stop[2] & 
                       e$chr == tads.sox5$chr[2])), ] 

e.sox9 <- e[which((e$bp.mid >= tads.sox9$bp.start[1] & 
                     e$bp.mid <= tads.sox9$bp.stop[1] & 
                     e$chr == tads.sox9$chr[1]) |
                    (e$bp.start >= tads.sox9$bp.start[1] & 
                       e$bp.start <= tads.sox9$bp.stop[1] & 
                       e$chr == tads.sox9$chr[1]) |
                    (e$bp.stop >= tads.sox9$bp.start[1] & 
                       e$bp.stop <= tads.sox9$bp.stop[1] & 
                       e$chr == tads.sox9$chr[1]) |
                    (e$bp.start >= tads.sox9$bp.start[2] & 
                       e$bp.start <= tads.sox9$bp.stop[2] & 
                       e$chr == tads.sox9$chr[2]) |
                    (e$bp.stop >= tads.sox9$bp.start[2] & 
                       e$bp.stop <= tads.sox9$bp.stop[2] & 
                       e$chr == tads.sox9$chr[2]) |
                    (e$bp.mid >= tads.sox9$bp.start[2] & 
                       e$bp.mid <= tads.sox9$bp.stop[2] & 
                       e$chr == tads.sox9$chr[2])), ] 

e.col10 <- e[which((e$bp.mid >= tads.col10$bp.start[1] & 
                     e$bp.mid <= tads.col10$bp.stop[1] & 
                     e$chr == tads.col10$chr[1]) |
                    (e$bp.start >= tads.col10$bp.start[1] & 
                       e$bp.start <= tads.col10$bp.stop[1] & 
                       e$chr == tads.col10$chr[1]) |
                    (e$bp.stop >= tads.col10$bp.start[1] & 
                       e$bp.stop <= tads.col10$bp.stop[1] & 
                       e$chr == tads.col10$chr[1]) |
                    (e$bp.start >= tads.col10$bp.start[2] & 
                       e$bp.start <= tads.col10$bp.stop[2] & 
                       e$chr == tads.col10$chr[2]) |
                    (e$bp.stop >= tads.col10$bp.start[2] & 
                       e$bp.stop <= tads.col10$bp.stop[2] & 
                       e$chr == tads.col10$chr[2]) |
                    (e$bp.mid >= tads.col10$bp.start[2] & 
                       e$bp.mid <= tads.col10$bp.stop[2] & 
                       e$chr == tads.col10$chr[2])), ]

rm(e)






cell.type.col <- c("orange", "#4682b4")

g1 <- ggplot()+
  geom_errorbarh(data = tads.sox5, aes(y = id, xmin = bp.start/10^6, xmax = bp.stop/10^6),
                 color = "darkgray", height = 0.25)+
  geom_errorbarh(data = genes.sox5, aes(y = "Genes", xmin = bp.start/10^6, xmax = bp.stop/10^6),
                 color = "black", height = 0.25)+
  geom_point(data = genes.sox5, aes(y = "Genes", x = tss/10^6), color = "red", size = 1)+
  geom_text_repel(data = genes.sox5, aes(y = "Genes", x = tss/10^6, label = Gene.name), 
                  color = "darkred", size = 3, min.segment.length = 0, nudge_y = 0.15)+
  geom_point(data = vista.sox5, aes(y = "Vista", x = bp.mid/10^6), 
             color = "black", shape = 1, size = 1)+
  geom_errorbarh(data = SE.sox5[SE.sox5$cell.type == "PC", ], aes(y = "SEs (PC)", xmin = bp.start/10^6, xmax = bp.stop/10^6), 
                 col = cell.type.col[1], height = 0.25)+
  geom_jitter(data = e.sox5[e.sox5$cell.type == "PC", ], aes(y = "H3K27ac (PC)", x = bp.mid/10^6), 
              color = cell.type.col[1], shape = 1, size = 1, height = 0.2, width = 0)+
  # geom_errorbarh(data = SE.sox5[SE.sox5$cell.type == "HC", ],aes(y = "SEs (HC)", xmin = bp.start/10^6, xmax = bp.stop/10^6),
  #                col = cell.type.col[2], height = 0.25)+
  geom_jitter(data = e.sox5[e.sox5$cell.type == "HC", ], aes(y = "H3K27ac (HC)", x = bp.mid/10^6), 
              color = cell.type.col[2], shape = 1, size = 1, height = 0.2, width = 0)+
  theme_bw(base_size = 12)+
  theme(legend.position = "top")+
  scale_color_manual(name = "Cell type", values = cell.type.col)+
  ylab(label = '')+
  xlab(label = 'Mb')+
  scale_y_discrete(limits = c("Vista", "H3K27ac (HC)",  "SEs (HC)", 
                              "H3K27ac (PC)", "SEs (PC)", "Genes", 
                              "TAD: HindIII", "TAD: NcoI"))+
  ggtitle(label = "Sox5 (chr 6)")




g2 <- ggplot()+
  geom_errorbarh(data = tads.sox9, aes(y = id, xmin = bp.start/10^6, xmax = bp.stop/10^6),
                 color = "darkgray", height = 0.25)+
  geom_errorbarh(data = genes.sox9, aes(y = "Genes", xmin = bp.start/10^6, xmax = bp.stop/10^6),
                 color = "black", height = 0.25)+
  geom_point(data = genes.sox9, aes(y = "Genes", x = tss/10^6), color = "red", size = 1)+
  geom_text_repel(data = genes.sox9, aes(y = "Genes", x = tss/10^6, label = Gene.name), 
                  color = "darkred", size = 3, min.segment.length = 0, nudge_y = 0.15)+
  geom_point(data = vista.sox9, aes(y = "Vista", x = bp.mid/10^6), 
             color = "black", shape = 1, size = 1)+
  geom_errorbarh(data = SE.sox9[SE.sox9$cell.type == "PC", ], aes(y = "SEs (PC)", xmin = bp.start/10^6, xmax = bp.stop/10^6), 
                 col = cell.type.col[1], height = 0.25)+
  geom_jitter(data = e.sox9[e.sox9$cell.type == "PC", ], aes(y = "H3K27ac (PC)", x = bp.mid/10^6), 
              color = cell.type.col[1], shape = 1, size = 1, height = 0.2, width = 0)+
  # geom_errorbarh(data = SE.sox9[SE.sox9$cell.type == "HC", ],aes(y = "SEs (HC)", xmin = bp.start/10^6, xmax = bp.stop/10^6),
  #                col = cell.type.col[2], height = 0.25)+
  geom_jitter(data = e.sox9[e.sox9$cell.type == "HC", ], aes(y = "H3K27ac (HC)", x = bp.mid/10^6), 
              color = cell.type.col[2], shape = 1, size = 1, height = 0.2, width = 0)+
  theme_bw(base_size = 12)+
  theme(legend.position = "top")+
  scale_color_manual(name = "Cell type", values = cell.type.col)+
  ylab(label = '')+
  xlab(label = 'Mb')+
  scale_y_discrete(limits = c("Vista", "H3K27ac (HC)",  "SEs (HC)", 
                              "H3K27ac (PC)", "SEs (PC)", "Genes", 
                              "TAD: HindIII", "TAD: NcoI"))+
  ggtitle(label = "Sox9 (chr 11)")



g <- gridExtra::arrangeGrob(g1, g2)
plot(g)
ggsave(filename = paste0("output/enhancers/Sox5_Sox9_plots.tiff"),
       plot = g, device = "tiff", width = 6, height = 6, dpi = 300)



g3 <- ggplot()+
  geom_errorbarh(data = tads.col10, aes(y = id, xmin = bp.start/10^6, xmax = bp.stop/10^6),
                 color = "darkgray", height = 0.25)+
  geom_errorbarh(data = genes.col10, aes(y = "Genes", xmin = bp.start/10^6, xmax = bp.stop/10^6),
                 color = "black", height = 0.25)+
  geom_point(data = genes.col10, aes(y = "Genes", x = tss/10^6), color = "red", size = 1)+
  geom_text_repel(data = genes.col10, aes(y = "Genes", x = tss/10^6, label = Gene.name), 
                  color = "darkred", size = 3, min.segment.length = 0, force = 3)+
  geom_point(data = vista.col10, aes(y = "Vista", x = bp.mid/10^6), 
             color = "black", shape = 1, size = 1)+
  # geom_errorbarh(data = SE.col10[SE.col10$cell.type == "PC", ], aes(y = "SEs (PC)", xmin = bp.start/10^6, xmax = bp.stop/10^6), 
  #                col = cell.type.col[1], height = 0.25)+
  geom_jitter(data = e.col10[e.col10$cell.type == "PC", ], aes(y = "H3K27ac (PC)", x = bp.mid/10^6), 
              color = cell.type.col[1], shape = 1, size = 1, height = 0.2, width = 0)+
  geom_errorbarh(data = SE.col10[SE.col10$cell.type == "HC", ],aes(y = "SEs (HC)", xmin = bp.start/10^6, xmax = bp.stop/10^6),
                 col = cell.type.col[2], height = 0.25)+
  geom_jitter(data = e.col10[e.col10$cell.type == "HC", ], aes(y = "H3K27ac (HC)", x = bp.mid/10^6), 
              color = cell.type.col[2], shape = 1, size = 1, height = 0.2, width = 0)+
  theme_bw(base_size = 12)+
  theme(legend.position = "top")+
  scale_color_manual(name = "Cell type", values = cell.type.col)+
  ylab(label = '')+
  xlab(label = 'Mb')+
  scale_y_discrete(limits = c("Vista", "H3K27ac (HC)",  "SEs (HC)", 
                              "H3K27ac (PC)", "SEs (PC)", "Genes", 
                              "TAD: HindIII", "TAD: NcoI"))+
  ggtitle(label = "Col10a1 (chr 10)")
ggsave(filename = paste0("output/enhancers/Col10a1_plots.tiff"),
       plot = g3, device = "tiff", width = 6, height = 3, dpi = 300)



g <- gridExtra::arrangeGrob(g1, g2, g3)
plot(g)
ggsave(filename = paste0("output/enhancers/Sox5_Sox9_Col10a1_plots.tiff"),
       plot = g, device = "tiff", width = 6, height = 9, dpi = 300)
