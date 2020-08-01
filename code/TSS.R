# Computes density of histone modification 
# peaks and ChromHMM states around TSS of 
# protein coding genes from Ensembl

require("ggplot2")
require("ggpubr")
require("IRanges")
require("parallel")
source("code/util/ChromHmmUtilities.R")
source("code/util/GeneralUtilities.R")

# fixed
cell.type.col <- c("orange", "#4682b4")


# Density of marks around TSS
getTssDensityMarks <- function() {
  
  absMin <- function(x) {
    return(min(abs(x)))
  }
  
  getAbsMin <- function(x) {
    x.dash <- numeric(length = length(x))
    for(i in 1:length(x)) {
      x.dash[i] <- absMin(x[i])
    }
    return(x[which.min(x.dash)[1]])
  }
  
  getTssDensityMarkSpecific <- function(i, marks) {
    mark <- marks[i]
    
    # data
    PC.peak.tss <- c()
    HC.peak.tss <- c()
    
    cat(mark, "\n")
    
    # PC
    file.PC <- list.files(path = "input/hidden_domain/PC/", 
                          pattern = mark, full.names = T)
    PC.peaks <- read.csv(file = file.PC, sep = "\t", 
                         as.is = TRUE, header = FALSE)
    colnames(PC.peaks) <- c("chr", "bp.start", "bp.stop")
    PC.peaks$bp.mid <- PC.peaks$bp.start+(PC.peaks$bp.stop-PC.peaks$bp.start)/2
    
    # HC
    file.HC <- list.files(path = "input/hidden_domain/HC/", 
                          pattern = mark, full.names = T)
    HC.peaks <- read.csv(file = file.HC, sep = "\t",
                         as.is = TRUE, header = FALSE)
    colnames(HC.peaks) <- c("chr", "bp.start", "bp.stop")
    HC.peaks$bp.mid <- HC.peaks$bp.start+(HC.peaks$bp.stop-HC.peaks$bp.start)/2
    
    # get unique chromosomes to study
    chrs <- unique(c(PC.peaks$chr, HC.peaks$chr))
    
    for(chr in chrs) {
      cat(chr, "\n")
      
      # TSS data
      tss.data <- getGenes(dir = "input/genes_biomart/genes.full.txt")
      tss.data <- tss.data[tss.data$Gene.type == "protein_coding", ]
      tss.data <- tss.data[tss.data$chr == chr, ]
      
      # PC
      temp <- PC.peaks[PC.peaks$chr == chr, ]
      if(nrow(temp) != 0) {
        temp$strand <- NA
        temp$min.tss.dist <- NA
        for(i in 1:nrow(temp)) {
          
          d <- temp$bp.mid[i]-tss.data$tss
          m <- which.min(abs(d))
          temp$min.tss.dist[i] <- d[m]
          temp$strand[i] <- tss.data$Strand[m]
          temp$mark <- mark
        }
        PC.peak.tss <- rbind(PC.peak.tss, temp)
      }
      
      # HC
      temp <- HC.peaks[HC.peaks$chr == chr, ]
      if(nrow(temp) != 0) {
        temp$strand <- NA
        temp$min.tss.dist <- NA
        for(i in 1:nrow(temp)) {
          
          d <- temp$bp.mid[i]-tss.data$tss
          m <- which.min(abs(d))
          temp$min.tss.dist[i] <- d[m]
          temp$strand[i] <- tss.data$Strand[m]
          temp$mark <- mark
        }
        HC.peak.tss <- rbind(HC.peak.tss, temp)
      }
    }
    rm(mark, d, i, temp, chr)
    
    
    # format correctly
    PC.peak.tss$cell.type <- "PC"
    HC.peak.tss$cell.type <- "HC"
    peak.tss <- rbind(PC.peak.tss, HC.peak.tss)
    peak.tss$min.tss.dist.corrected <- ifelse(
      test = peak.tss$strand == -1,
      yes = peak.tss$min.tss.dist*(-1),
      no = peak.tss$min.tss.dist*1)
    
    
    peak.tss$cell.type <- factor(x = peak.tss$cell.type, 
                                 levels = c("PC", "HC"))
    peak.tss$mark <- factor(x = peak.tss$mark, 
                            levels = getMarks())
    
    return (peak.tss)
  }
  
  
  marks <- getMarks()
  
  peak.tss <- mclapply(X = 1:length(marks),
                       FUN = getTssDensityMarkSpecific,
                       marks = marks,
                       mc.cores = 10)
  peak.tss <- do.call(rbind, peak.tss)
  
  return (peak.tss)
}



# Density of states around TSS
getTssDensityState <- function() {
  
  absMin <- function(x) {
    return(min(abs(x)))
  }
  
  getAbsMin <- function(x) {
    x.dash <- numeric(length = length(x))
    for(i in 1:length(x)) {
      x.dash[i] <- absMin(x[i])
    }
    return(x[which.min(x.dash)[1]])
  }
  
  # PC
  ePC <- get(load(file = "output/enhancers/ePC.RData"))
  ePC <- ePC[, c("chr", "bp.start", "bp.stop")]
  
  # HC
  eHC <- get(load(file = "output/enhancers/eHC.RData"))
  eHC <- eHC[, c("chr", "bp.start", "bp.stop")]
  
  
  ePC.tss <- c()
  eHC.tss <- c()
  for(chr in unique(c(ePC$chr, eHC$chr))) {
    cat(chr, "\n")
    
    # TSS data
    tss.data <- getGenes(dir = "input/genes_biomart/genes.full.txt")
    tss.data <- tss.data[tss.data$Gene.type == "protein_coding", ]
    tss.data <- tss.data[tss.data$chr == chr, ]
    
    # PC
    temp <- ePC[ePC$chr == chr, ]
    if(nrow(temp) != 0) {
      temp$min.tss.dist <- NA
      temp$strand <- NA
      for(i in 1:nrow(temp)) {
        
        # if bp.start:bp.stop around a TSS => TRUE => 0 distance, else => INF distance
        o <- (temp$bp.start[i] <= tss.data$tss & temp$bp.stop[i] >= tss.data$tss)
        o <- ifelse(test = o == TRUE, yes = 0, no = Inf)
        
        d <- cbind(temp$bp.start[i]-tss.data$tss,
                   temp$bp.stop[i]-tss.data$tss,
                   o)
        
        m <- which.min(apply(X = d, MARGIN = 1, FUN = absMin))
        temp$min.tss.dist[i] <- getAbsMin(x = d[m, ])
        temp$strand[i] <- tss.data$Strand[m]
      }
      ePC.tss <- rbind(ePC.tss, temp)
    }
    
    # HC
    temp <- eHC[eHC$chr == chr, ]
    if(nrow(temp) != 0) {
      temp$min.tss.dist <- NA
      temp$strand <- NA
      for(i in 1:nrow(temp)) {
        
        # if bp.start:bp.stop around a TSS => TRUE => 0 distance, else => INF distance
        o <- (temp$bp.start[i] <= tss.data$tss & temp$bp.stop[i] >= tss.data$tss)
        o <- ifelse(test = o == TRUE, yes = 0, no = Inf)
        
        d <- cbind(temp$bp.start[i]-tss.data$tss,
                   temp$bp.stop[i]-tss.data$tss,
                   o)
        
        m <- which.min(apply(X = d, MARGIN = 1, FUN = absMin))
        temp$min.tss.dist[i] <- getAbsMin(x = d[m, ])
        temp$strand[i] <- tss.data$Strand[m]
      }
      eHC.tss <- rbind(eHC.tss, temp)
    }
  }
  rm(temp, i, d, chr)
  
  ePC.tss$cell.type <- "PC"
  eHC.tss$cell.type <- "HC"
  peak.tss <- rbind(ePC.tss, eHC.tss)
  peak.tss$peak.count <- 1
  peak.tss$min.tss.dist.corrected <- ifelse(test = peak.tss$strand == -1,
                                            yes = peak.tss$min.tss.dist*(-1),
                                            no = peak.tss$min.tss.dist*1)
  rm(ePC.tss, eHC.tss, ePC, eHC)
  
  peak.tss$cell.type <- factor(x = peak.tss$cell.type, 
                               levels = c("PC", "HC"))
  
  return (peak.tss)
}



##### Density of peaks around TSS ######

peak.tss <- getTssDensityMarks()
save(peak.tss, file = "output/tss/mark.tss.RData")

g <- gridExtra::arrangeGrob(
  ncol = 2,
  ggplot(data = peak.tss[peak.tss$min.tss.dist.corrected > -5*10^3 & 
                           peak.tss$min.tss.dist.corrected < 5*10^3, ])+
    facet_wrap(facets = ~mark, ncol = 1, scale = "free_y")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_density(aes(min.tss.dist.corrected/10^3, col = cell.type))+
    theme_bw(base_size = 11)+
    theme(legend.position = "top")+
    ylab(label = "Density of peaks")+
    xlab(label = "Distance to TSS [Kb]")+
    scale_color_manual(name = "Cell type", 
                       values = cell.type.col),
  ggplot(data = peak.tss[peak.tss$min.tss.dist.corrected > -5*10^4 & 
                           peak.tss$min.tss.dist.corrected < 5*10^4, ])+
    facet_wrap(facets = ~mark, ncol = 1, scale = "free_y")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_density(aes(min.tss.dist.corrected/10^3, col = cell.type))+
    theme_bw(base_size = 11)+
    theme(legend.position = "top")+
    ylab(label = "Density of peaks")+
    xlab(label = "Distance to TSS [Kb]")+
    scale_color_manual(name = "Cell type", 
                       values = cell.type.col))
plot(g)

ggsave(plot = g, filename = "output/tss/Mark_tss.tiff",
       device = "tiff", width = 6, height = 8, dpi = 300)



# discretized density of peaks in ranges
ranges <- c(-10^5, -5*10^4, -10^4, -5*10^3, 5*10^3, 10^4, 5*10^4, 10^5)
range.keys <- c()
peak.tss$range <- NA

for(i in 1:(length(ranges)-1)) {
  key <- paste("[", paste(ranges[i]/10^3, ranges[i+1]/10^3, sep = ', '), "]", sep = '')
  peak.tss$range[which(peak.tss$min.tss.dist.corrected >= ranges[i] & 
                         peak.tss$min.tss.dist.corrected <= ranges[i+1])] <- key
  range.keys <- c(range.keys, key)
}
peak.tss$peak.count <- 1

summarized.peak.tss <- aggregate(peak.count~cell.type+range+mark, 
                                 data = peak.tss, FUN = sum)
summarized.peak.tss$range <- factor(x = summarized.peak.tss$range, 
                                    levels = range.keys)


g <- ggplot(data = summarized.peak.tss)+
  facet_wrap(facets = ~mark, ncol = 2)+
  geom_bar(aes(x = range, weight = peak.count, fill = cell.type),
           position = position_dodge(width = 0.5), width = 0.5)+
  theme_bw(base_size = 11)+
  theme(legend.position = "top")+
  scale_fill_manual(name = "Cell type", values = c("#4682b4", "orange"))+
  xlab(label = "Ranges [Kb]")+
  ylab(label = "Peak count")
g
ggsave(plot = g, filename = "output/tss/Mark_bins_tss.tiff",
       device = "tiff", width = 8.5, height = 5, dpi = 300)





##### Density of states around ChromHMM state 7 ######

peak.tss <- getTssDensityState()

g <- gridExtra::arrangeGrob(
  ncol = 2,
  ggplot(data = peak.tss[peak.tss$min.tss.dist.corrected > -5*10^3 & 
                           peak.tss$min.tss.dist.corrected < 5*10^3, ])+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_density(aes(min.tss.dist.corrected/10^3, col = cell.type))+
    theme_bw(base_size = 11)+
    theme(legend.position = "top")+
    ylab(label = "Density of peaks")+
    xlab(label = "Distance to TSS [Kb]")+
    scale_color_manual(name = "Cell type", 
                       values = cell.type.col),
  ggplot(data = peak.tss[peak.tss$min.tss.dist.corrected > -5*10^4 & 
                           peak.tss$min.tss.dist.corrected < 5*10^4, ])+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_density(aes(min.tss.dist.corrected/10^3, col = cell.type))+
    theme_bw(base_size = 11)+
    theme(legend.position = "top")+
    ylab(label = "Density of peaks")+
    xlab(label = "Distance to TSS [Kb]")+
    scale_color_manual(name = "Cell type", 
                       values = cell.type.col))

plot(g)

ggsave(plot = g, filename = "output/tss/State7_tss.tiff",
       device = "tiff", width = 6, height = 3, dpi = 300)



# discretized density of peaks in ranges
ranges <- c(-10^5, -5*10^4, -10^4, -5*10^3, 5*10^3, 10^4, 5*10^4, 10^5)
range.keys <- c()
peak.tss$range <- NA

for(i in 1:(length(ranges)-1)) {
  key <- paste("[", paste(ranges[i]/10^3, ranges[i+1]/10^3, sep = ', '), "]", sep = '')
  peak.tss$range[which(peak.tss$min.tss.dist.corrected >= ranges[i] & 
                         peak.tss$min.tss.dist.corrected <= ranges[i+1])] <- key
  range.keys <- c(range.keys, key)
}
peak.tss$peak.count <- 1

summarized.peak.tss <- aggregate(peak.count~cell.type+range, 
                                 data = peak.tss, FUN = sum)
summarized.peak.tss$range <- factor(x = summarized.peak.tss$range, 
                                    levels = range.keys)


g <- ggplot(data = summarized.peak.tss)+
  geom_bar(aes(x = range, weight = peak.count, fill = cell.type),
           position = position_dodge(width = 0.5), width = 0.5)+
  theme_bw(base_size = 11)+
  theme(legend.position = "top")+
  scale_fill_manual(name = "Cell type", values = c("#4682b4", "orange"))+
  xlab(label = "Ranges [Kb]")+
  ylab(label = "Peak count")
g
ggsave(plot = g, filename = "output/tss/State7_bins_tss.tiff",
       device = "tiff", width = 5, height = 4, dpi = 300)
