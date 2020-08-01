# Input: HiddenDomain
# Output: input for ChromHMM

source("code/util/ChromHmmUtilities.R")
source("code/util/GeneralUtilities.R")
options(scipen=999)

# split genome into bins
chr.len <- read.csv(file = "input/ChromHMM/CHROMSIZES/mm10.txt", 
                    sep = "\t", header = F, as.is = T)
colnames(chr.len) <- c("chr", "length")
chr.len <- chr.len[chr.len$chr %in% getChrs(), ]


# const
bin.size <- 200


chrs <- getChrs()
header <- c()
for(chr in chrs) {
  cat(chr, "\n")
  j <- which(chr.len$chr == chr)
  temp <- data.frame(chr = chr.len$chr[j], stringsAsFactors = F,
                     bp.start = seq(from = 0, to = chr.len$length[j], by = bin.size))
  temp$bp.stop <- temp$bp.start + bin.size
  header <- rbind(header, temp)
}
rm(chr, j, chrs, temp)
save(header, file = "output/header.RData")


write.table(x = header, file = "output/hidden_domain_bins/header.bed", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = F)




# b. bamming
files <- c(list.files(path = "input/hidden_domain/PC/", pattern = ".bed", 
                      all.files = T, recursive = F, full.names = T),
           list.files(path = "input/hidden_domain/HC/", pattern = ".bed", 
                      all.files = T, recursive = F, full.names = T))
short.files <- c(list.files(path = "input/hidden_domain/PC/", pattern = ".bed", 
                            all.files = T, recursive = F, full.names = F),
                 list.files(path = "input/hidden_domain/HC/", pattern = ".bed", 
                            all.files = T, recursive = F, full.names = F))

short.files <- do.call(rbind, strsplit(x = short.files, split = "_Intersect"))[, 1]
for(i in 1:length(files)) {
  system(command = paste("bedtools intersect -a output/hidden_domain_bins/header.bed -b ", 
                         files[i], " -u > output/hidden_domain_bins/", short.files[i], sep = ''))
  cat(i, "\n")
}
rm(i)




# c. merge data
header <- read.csv(file = "output/hidden_domain_bins/header.bed", sep = "\t", header = F, as.is = T)
header$key <- paste(header[, 1], header[, 2], sep = '.')
files <- list.files(path = "output/hidden_domain_bins/", pattern = "HC|PC", full.names = T)
short.files <- list.files(path = "output/hidden_domain_bins/", pattern = "HC|PC", full.names = F)

for(i in 1:length(files)) {
  o <- read.csv(file = files[i], stringsAsFactors = F, sep = "\t", header = F)
  o$V3 <- NULL
  colnames(o) <- c("chr", "bp")
  o$key <- paste(o$chr, o$bp, sep = '.')
  hits <- which(header$key %in% o$key)
  
  cell.state <- ifelse(test = regexpr(pattern = "HC", text = short.files[i], 
                                      ignore.case = F) == -1, yes = "PC", no = "HC")
  mark.state <- gsub(pattern = "PC-|HC-", replacement = '', x = short.files[i])
  header[, paste(cell.state, mark.state, sep = '.')] <- 0
  header[hits, paste(cell.state, mark.state, sep = '.')] <- 1
}
header$key <- NULL
hm.bins <- header
colnames(hm.bins)[1:3] <- c("chr", "bp.start", "bp.stop")
save(hm.bins, file = "output/hidden_domain_bins/hm.bins.RData")
rm(i, o, hits, cell.state, mark.state)





# d. Prepare hidden domain data for ChromHmm
del.chr <- paste0(paste(paste0("S", 1:100), collapse = "|"), 
                  paste(paste0("ID", 1:100), collapse = "|"),
                  "|\\.bed|-|Intersect|PC\\.|HC\\.", collapse = '')

hm.bins <- get(load(file = "output/hidden_domain_bins/hm.bins.RData"))

chromosomes <- getChrs()
chromosomes <- chromosomes[chromosomes != "chrMt"]
for(chr in chromosomes) {
  temp <- hm.bins[hm.bins$chr == chr, ]
  
  # PC for chromhmm
  pc.i <- which(regexpr(pattern = "PC\\.", text = colnames(temp), ignore.case = F) != -1)
  pc <- temp[, pc.i]
  colnames(pc) <- gsub(pattern = del.chr, replacement = '', x = colnames(pc))
  pc <- rbind(colnames(pc), pc)
  pc <- rbind(c("PC", chr, rep(x = "", times = ncol(pc))), pc)
  write.table(x = pc, append = F, quote = F,  row.names = F, col.names = F, sep = "\t", 
              file = paste("output/chromhmm/PC_", chr, "_binary.csv", sep = ''))
  
  
  # HC for chromhmm
  hc.i <- which(regexpr(pattern = "HC\\.", text = colnames(temp), ignore.case = F) != -1)
  hc <- temp[, hc.i]
  colnames(hc) <- gsub(pattern = del.chr, replacement = '', x = colnames(hc))
  hc <- rbind(colnames(hc), hc)
  hc <- rbind(c("HC", chr, rep(x = "", times = ncol(hc))), hc)
  write.table(x = hc, append = F, quote = F,  row.names = F, col.names = F, sep = "\t", 
              file = paste("output/chromhmm/HC_", chr, "_binary.csv", sep = ''))
  
  cat(chr, "\n")
}
