# Correlation analysis: ChromHMM state prevalence ~ Gene expression



# requirements
source("code/util/GeneralUtilities.R")
source("code/util/ChromHmmUtilities.R")
require("boot")
require("reshape2")
require("ranger")
require("ggplot2")
require("parallel")


#### Prepare data for analysis ####

# gene-centric data
out <- get(load(file = "output/gene_centric/gcsc3.states.RData"))
out$state <- as.character(out$state)


l <- table(out$gene.id, out$state, out$region, out$cell.type)
l <- data.frame(l, stringsAsFactors = F)
colnames(l) <- c("gene.id", "state", "region", "cell.type", "count.cor")
l$gene.id <- as.character(l$gene.id)
l$state <- as.character(l$state)
l$region <- as.character(l$region)
l$cell.type <- as.character(l$cell.type)

out <- merge(x = l, y = out, all.x = T,
           by = c("gene.id", "state", 
                  "region", "cell.type"))
out$strand <- NULL
out$chr <- NULL
out$count <- out$count.cor
out$count.cor <- NULL
write.table(x = out, file = "output/tpm_states/out.tsv", append = F,
            quote = F, sep = "\t", row.names = F, col.names = T)
rm(l, o)





#### Analysis ####
# requirements
source("code/util/GeneralUtilities.R")
source("code/util/ChromHmmUtilities.R")
require("boot")
require("reshape2")
require("ranger")
require("ggplot2")
require("parallel")


# fixed
cell.type.col <- c("orange", "#4682b4")


# Description:
getPrevalence <- function(i, ranks, de.tpm, out, 
                          hdi.level, boot.N) {
  
  # boot mean
  boot.mean <- function(data, index) {
    return(mean(data[index]))
  }
  
  
  stats <- c()
  for(region in c("enhancer", "tss", "gene")) {
    for(state in 1:15) {
      
      # PC data
      temp.tpm <- de.tpm[de.tpm$rank.P >= (ranks[i]+1) & 
                           de.tpm$rank.P <= ranks[i+1], ]
      
      temp.gcsc <- out[out$region == region & 
                         out$cell.type == "PC" & 
                         out$gene.id %in% temp.tpm$ens_gene &
                         out$state == state, ]
      
      b <- boot(data = temp.tpm$mean.P, statistic = boot.mean, R = boot.N)
      mu <- mean(b$t[, 1])
      hdi <- getHdi(vec = b$t[, 1], hdi.level = hdi.level)
      
      b.p <- boot(data = ifelse(test = temp.gcsc$count == 0, yes = 0, no = 1), 
                  statistic = boot.mean, R = boot.N)
      mu.p <- mean(b.p$t[, 1])
      hdi.p <- getHdi(vec = b.p$t[, 1], hdi.level = hdi.level)
      
      row <- data.frame(rank = i, rank.min = ranks[i]+1, rank.max = ranks[i+1],
                        region = region, state = state, cell.type = "PC",
                        tpm.mean = mu, tpm.hdi.95.L = hdi[1], tpm.hdi.95.H = hdi[2],
                        prevalence.mean = mu.p, prevalence.hdi.95.L = hdi.p[1], 
                        prevalence.hdi.95.H = hdi.p[2], stringsAsFactors = FALSE)
      stats <- rbind(stats, row)
      
      
      # HC data
      temp.tpm <- de.tpm[de.tpm$rank.H >= (ranks[i]+1) & 
                           de.tpm$rank.H <= ranks[i+1], ]
      
      temp.gcsc <- out[out$region == region & 
                         out$cell.type == "HC" & 
                         out$gene.id %in% temp.tpm$ens_gene &
                         out$state == state, ]
      
      b <- boot(data = temp.tpm$mean.H, statistic = boot.mean, R = boot.N)
      mu <- mean(b$t[, 1])
      hdi <- getHdi(vec = b$t[, 1], hdi.level = hdi.level)
      
      b.p <- boot(data = ifelse(test = temp.gcsc$count == 0, yes = 0, no = 1), 
                  statistic = boot.mean, R = boot.N)
      mu.p <- mean(b.p$t[, 1])
      hdi.p <- getHdi(vec = b.p$t[, 1], hdi.level = hdi.level)
      
      row <- data.frame(rank = i, rank.min = ranks[i]+1, rank.max = ranks[i+1],
                        region = region, state = state, cell.type = "HC",
                        tpm.mean = mu, tpm.hdi.95.L = hdi[1], tpm.hdi.95.H = hdi[2],
                        prevalence.mean = mu.p, prevalence.hdi.95.L = hdi.p[1], 
                        prevalence.hdi.95.H = hdi.p[2], stringsAsFactors = FALSE)
      stats <- rbind(stats, row)
    }
  }
  return (stats)
}


# Description:
getPrevalenceZero <- function(t, hdi.level, boot.N, states) {
  
  # boot mean
  boot.mean <- function(data, index) {
    return(mean(data[index]))
  }
  
  regions <- unique(t$region)
  cell.types <- unique(t$cell.type)
  
  stats <- c()
  for(r in regions) {
    for(s in states) {
      for(c in cell.types) {
        
        d <- t[which(t$region == r & t$state == s & t$cell.type == c), ]
        b.p <- boot(data = d$count, statistic = boot.mean, R = boot.N)
        mu.p <- mean(b.p$t[, 1])
        hdi.p <- getHdi(vec = b.p$t[, 1], hdi.level = hdi.level)
        
        row <- data.frame(rank = 0, rank.min = 0, rank.max = 0,
                          region = r, state = s, cell.type = c,
                          tpm.mean = 0, tpm.hdi.95.L = 0, tpm.hdi.95.H = 0,
                          prevalence.mean = mu.p, prevalence.hdi.95.L = hdi.p[1], 
                          prevalence.hdi.95.H = hdi.p[2], stringsAsFactors = FALSE)
        stats <- rbind(stats, row)
      }
    }
  }
  
  return (stats)
}


# Description:
# TRUE if all(x == 0)
allZero <- function(x) {
  return(all(x==0))
}



out.data <- read.csv(file = "output/tpm_states/out.tsv",
                     header = TRUE, sep = "\t", as.is = TRUE)



# rna-seq data
de.tpm <- read.csv(file = "output/tpms_gene.csv", sep = ",",  
                   header = T, as.is = T)
de.tpm <- de.tpm[which(apply(X = de.tpm[, c("P1", "P2", "H3", "H4"), ], 
                             MARGIN = 1, FUN = allZero) == F), ]
de.tpm$mean.P <- apply(X = de.tpm[, c("P1", "P2")], MARGIN = 1, FUN = mean)
de.tpm$sd.P <- apply(X = de.tpm[, c("P1", "P2")], MARGIN = 1, FUN = sd)
de.tpm$mean.H <- apply(X = de.tpm[, c("H3", "H4")], MARGIN = 1, FUN = mean)
de.tpm$sd.H <- apply(X = de.tpm[, c("H3", "H4")], MARGIN = 1, FUN = sd)


out <- out.data[out.data$gene.id %in% de.tpm$ens_gene, ]
de.tpm <- de.tpm[de.tpm$ens_gene %in% out$gene.id, ]


# rank P and H by mean TPM
de.tpm <- de.tpm[order(de.tpm$mean.P, decreasing = T), ]
de.tpm$rank.P <- 1:nrow(de.tpm)
de.tpm <- de.tpm[order(de.tpm$mean.H, decreasing = T), ]
de.tpm$rank.H <- 1:nrow(de.tpm)


ranks <- seq(from = 0, to = nrow(de.tpm), by = 500)
ranks[length(ranks)] <- nrow(de.tpm)
stats <- mclapply(X = 1:(length(ranks)-1), 
                  FUN = getPrevalence, 
                  ranks = ranks, 
                  de.tpm, out, 
                  hdi.level = 0.95, 
                  boot.N = 10^4,
                  mc.cores = 20)
stats <- do.call(rbind, stats)
save(stats, file = "output/tpm_states/stats.15.RData")



# plot
stats <- get(load("~/Chrom/output/tpm_states/stats.15.RData"))
stats$cell.type <- factor(x = stats$cell.type, levels = c("PC", "HC"))

g <- ggplot(data = stats)+
  facet_grid(facets = state~region)+
  geom_point(aes(x = tpm.mean, y = prevalence.mean, col = cell.type),
             size = 1, shape = 21, fill = NA)+
  geom_errorbar(aes(x = tpm.mean, ymin = prevalence.hdi.95.L,
                    ymax = prevalence.hdi.95.H, col = cell.type),
                width = 0.03)+
  theme_bw(base_size = 12)+
  scale_x_continuous(trans = "log10",
                     breaks = c(10^seq(from=-4, to = 4, by = 2)))+
  annotation_logticks(base = 10, sides = "b")+
  scale_color_manual(name = "Cell state", values = cell.type.col)+
  scale_shape_manual(name = "Cell state", values = c(21, 22))+
  theme(legend.position = "top")+
  ylab("Mean state prevalence")+
  xlab("Mean TPM")+
  ylim(c(0, 1))
g
ggsave(plot = g, file = "output/tpm_states/scaling_updated.tiff",
       device = "tiff", width = 9, height = 16, dpi = 300)
ggsave(plot = g, file = "output/tpm_states/scaling_updated.pdf",
       device = "pdf", width = 9, height = 16, dpi = 300)



