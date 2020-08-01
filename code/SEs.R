# Super enhancers 


# required code/libs
source("code/util/GeneralUtilities.R")
source("code/util/ChromHmmUtilities.R")
source("code/util/TadStats.R")
source("code/util/EnhancerUtillities.R")
require(reshape)
require(parallel)
require(IRanges)



# Description:
# Get distances between enhancers
getEDistances <- function(d) {
  
  getClust <- function(x) {
    x <- x[order(x$bp.start, decreasing = F), ]
    d <- 0
    for(i in 2:nrow(x)) {
      d <- c(d, (abs(x$bp.start[i]-x$bp.stop[i-1])/10^3))
    }
    
    x$dist.kb <- d
    return (x)
  }
  
  out <- c()
  chrs <- unique(d$chr)
  for(chr in chrs) {
    out <- rbind(out, getClust(x = d[d$chr %in% chr, ]))
  }
  
  return (out)
}



# Description:
# Get super-enhancer (stitch)
getSeMax <- function(d, gap.kb.thr) {
  
  getClust <- function(x, gap.kb.thr) {
    x <- x[order(x$bp.start, decreasing = F), ]
    x$SE_ID <- NA
    x$SE_ID[1] <- 1
    d.cont <- T
    for(i in 2:nrow(x)) {
      d.cont <- ((abs(x$bp.start[i]-x$bp.stop[i-1])/10^3) <= gap.kb.thr)
      if(d.cont) {
        x$SE_ID[i] <- x$SE_ID[i-1]
      }
      else {
        x$SE_ID[i] <- x$SE_ID[i-1]+1
      }
    }
    
    return (x)
  }
  
  out <- c()
  chrs <- unique(d$chr)
  for(chr in chrs) {
    x <- getClust(x = d[d$chr %in% chr, ], gap.kb.thr = gap.kb.thr)
    out <- rbind(out, x)
  }
  return (out)
}



# load peaks data
PC <- read.csv(file = "input/hidden_domain/PC/PC-H3K27ac-S3-S5-S8-Intersect-ID35-ID39-ID47.bed",
               header = FALSE, sep = "\t", as.is = TRUE)
colnames(PC) <- c("chr", "bp.start", "bp.stop")
PC$bp.mid <- PC$bp.start + (PC$bp.stop - PC$bp.start)/2
PC$cell.type <- "PC"
PC$overlap.TSS <- 0

HC <- read.csv(file = "input/hidden_domain/HC/HC-H3K27ac-S5-S7-Intersect-ID19-ID25.bed",
               header = FALSE, sep = "\t", as.is = TRUE)
colnames(HC) <- c("chr", "bp.start", "bp.stop")
HC$bp.mid <- HC$bp.start + (HC$bp.stop - HC$bp.start)/2
HC$cell.type <- "HC"
HC$overlap.TSS <- 0



# low/high limit of TSS +/- 2.5kb
gene.data <- getGenes(dir = "input/genes_biomart/genes.full.txt")
gene.data <- gene.data[gene.data$Gene.type == "protein_coding", ]
gene.data$tss.L = gene.data$tss-2500
gene.data$tss.H = gene.data$tss+2500



# determine peaks that overlap TSS
ePC <- c()
eHC <- c()
chrs <- unique(c(PC$chr, HC$chr))
for(chr in chrs) {
  
  # PC
  temp.query <- PC[PC$chr == chr, ]
  temp.subject <- gene.data[gene.data$chr == chr, ]
  
  q <- IRanges(start = temp.query$bp.start, end = temp.query$bp.stop)
  s <- IRanges(start = temp.subject$tss.L, end = temp.subject$tss.H)
  o <- findOverlaps(query = q, subject = s, type = "any")
  o <- unique(o@from)
  temp.query$overlap.TSS[o] <- 1
  ePC <- rbind(ePC, temp.query)
  rm(o, q, s, temp.query, temp.subject)
  
  
  # HC
  temp.query <- HC[HC$chr == chr, ]
  temp.subject <- gene.data[gene.data$chr == chr, ]
  
  q <- IRanges(start = temp.query$bp.start, end = temp.query$bp.stop)
  s <- IRanges(start = temp.subject$tss.L, end = temp.subject$tss.H)
  o <- findOverlaps(query = q, subject = s, type = "any")
  o <- unique(o@from)
  temp.query$overlap.TSS[o] <- 1
  eHC <- rbind(eHC, temp.query)
  rm(o, q, s, temp.query, temp.subject)
}
rm(chr, chrs, PC, HC)
save(ePC, file = "output/enhancers/ePC.RData")
save(eHC, file = "output/enhancers/eHC.RData")



# remove peaks overlapping TSS +/- 2.5
sePC <- ePC[ePC$overlap.TSS == 0, ]
seHC <- eHC[eHC$overlap.TSS == 0, ]

# Tables: peaks non-TSS overlapping
sePC <- ePC[ePC$overlap.TSS == 0, ]
write.table(x = sePC[, 1:3], file = "output/enhancers/H3K27ac_PC_nonTSS.bam", 
            sep = "\t", append = F, quote = F, row.names = F, col.names = F)
seHC <- eHC[eHC$overlap.TSS == 0, ]
write.table(x = seHC[, 1:3], file = "output/enhancers/H3K27ac_HC_nonTSS.bam", 
            sep = "\t", append = F, quote = F, row.names = F, col.names = F)


# Distances and super-enhancer clustering 
sePC <- getEDistances(d = sePC)
seHC <- getEDistances(d = seHC)




# Set SE_IDs
sePC <- getSeMax(d = sePC, gap.kb.thr = 12.5)
sePC$SE_ID <- paste0("SE_", sePC$chr, "_", sePC$SE_ID)
seHC <- getSeMax(d = seHC, gap.kb.thr = 12.5)
seHC$SE_ID <- paste0("SE_", seHC$chr, "_", seHC$SE_ID)


# Aggregate SEs
SE <- rbind(sePC, seHC)
SE$count <- 1
rm(sePC, seHC)


SEmin <- aggregate(formula = bp.start~chr+cell.type+SE_ID, data = SE, FUN = min)
SEmax <- aggregate(formula = bp.stop~chr+cell.type+SE_ID, data = SE, FUN = max)
SEcount <- aggregate(formula = count~chr+cell.type+SE_ID, data = SE, FUN = sum)
SE.data <- merge(x = SEmax, y = SEmin, by = c("chr", "cell.type", "SE_ID"))
SE.data <- merge(x = SE.data, y = SEcount, by = c("chr", "cell.type", "SE_ID"))
rm(SEmin, SEmax, SEcount, SE)


SE.data <- SE.data[order(SE.data$count, decreasing = TRUE), ]
save(SE.data, file = "output/enhancers/SE.RData")

SE.data <- SE.data[which(SE.data$count >= 5), ]
# Tables: SEs
x <- SE.data[SE.data$cell.type == "PC", ]
write.table(x = x[, c("chr", "bp.start", "bp.stop")], 
            file = "output/enhancers/SE_PC.bam", 
            sep = "\t", append = F, quote = F, 
            row.names = F, col.names = F)

y <- SE.data[SE.data$cell.type == "HC", ]
write.table(x = y[, c("chr", "bp.start", "bp.stop")], 
            file = "output/enhancers/SE_HC.bam", 
            sep = "\t", append = F, quote = F, 
            row.names = F, col.names = F)

