



getMarkEnrichmentDistribution <- function(chr) { 
  # Col2s
  col2.files <- list.files(path = "output/deseq2/", pattern = paste(chr, "\\.Col2\\.", sep = ''))
  col2.names <- gsub(pattern = paste(chr, "\\.Col2\\.", sep = ''), replacement = '', x = col2.files)
  col2.names <- gsub(pattern = "\\.RData", replacement = '', x = col2.names)
  
  col2.lfc <- c()
  for(i in 1:length(col2.files)) {
    dseq.out <- get(load(paste("output/deseq2/", col2.files[i], sep = '')))
    temp <- data.frame(lfc = dseq.out$log2FoldChange, 
                       base.mean = dseq.out$baseMean, 
                       p.value = dseq.out$pvalue,
                       fdr = dseq.out$padj,
                       stat = dseq.out$stat,
                       mark = col2.names[i], 
                       state = "col2")
    temp$bit <- ifelse(test = temp$lfc >= 0, yes = 1, no = 0)
    col2.lfc <- rbind(col2.lfc, temp)
  }
  
  
  # Col10s
  col10.files <- list.files(path = "output/deseq2/", pattern = paste(chr, "\\.Col10\\.", sep = ''))
  col10.names <- gsub(pattern = paste(chr, "\\.Col10\\.", sep = ''), replacement = '', x = col10.files)
  col10.names <- gsub(pattern = "\\.RData", replacement = '', x = col10.names)
  
  col10.lfc <- c()
  for(i in 1:length(col10.files)) {
    dseq.out <- get(load(paste("output/deseq2/", col10.files[i], sep = '')))
    temp <- data.frame(lfc = dseq.out$log2FoldChange,
                       base.mean = dseq.out$baseMean,
                       p.value = dseq.out$pvalue,
                       fdr = dseq.out$padj,
                       mark = col10.names[i],
                       stat = dseq.out$stat,
                       state = "col10")
    temp$bit <- ifelse(test = temp$lfc >= 0, yes = 1, no = 0)
    col10.lfc <- rbind(col10.lfc, temp)
  }
  
  lfc <- rbind(col2.lfc, col10.lfc)
  return(lfc)
}





getEmissionsAndTransitions <- function(nr.states, 
                                       chrom.dir) {
  
  parseModelFile <- function(file) {
    lines <- readLines(con = file)
    emission.lines <- lines[which(regexpr(pattern = "emissionprobs", text = lines) != -1)]
    transition.lines <- lines[which(regexpr(pattern = "transitionprobs", text = lines) != -1)]
    
    rows <- c()
    for(i in 1:length(emission.lines)) {
      row <- unlist(strsplit(x = emission.lines[i], split = "\t"))
      row <- row[c(2, 4, 5, 6)]
      rows <- rbind(rows, row)
    }
    emissions <- data.frame(state = as.numeric(rows[, 1]), 
                            mark = rows[, 2], 
                            emission.bit = as.numeric(rows[, 3]),
                            emission = as.numeric(rows[, 4]))
    emissions <- emissions[emissions$emission.bit == 1, ]
    
    rows <- c()
    for(i in 1:length(transition.lines)) {
      row <- unlist(strsplit(x = transition.lines[i], split = "\t"))
      rows <- rbind(rows, row[2:4])
    }
    transitions <- data.frame(state.from = as.numeric(rows[, 1]),
                              state.to = as.numeric(rows[, 2]),
                              transition = as.numeric(rows[, 3]))
    
    emissions.m <- acast(data = emissions, value.var = "emission", 
                         formula = state ~ mark)
    transitions.m <- acast(data = transitions, value.var = "transition", 
                           formula = state.from ~ state.to)
    
    return(list(emissions = emissions,
                transitions = transitions, 
                emissions.m = emissions.m, 
                transitions.m = transitions.m))
  }
  
  
  parseOverlapFile <- function(file) {
    overlap <- read.csv(file = file, header = T, sep = "\t", stringsAsFactors = F)
    
    # overlap filter
    rownames(overlap) <- overlap$state..Emission.order.
    overlap$state..Emission.order. <- NULL
    overlap <- overlap
    
    # data
    overlap.data <- c()
    for(i in 1:nrow(overlap)) {
      for(j in 1:ncol(overlap)) {
        row <- c(rownames(overlap)[i], colnames(overlap)[j], overlap[i, j])
        overlap.data <- rbind(overlap.data, row)
      }  
    }
    
    overlap.data <- data.frame(state = overlap.data[, 1], 
                               region = overlap.data[, 2], 
                               overlap = as.numeric(gsub(pattern = '\\,', replacement = '.', x = overlap.data[, 3])), 
                               stringsAsFactors = F)
    
    return (list(overlap = overlap.data, overlap.m = overlap))
  }
  
  
  file.model <- paste(chrom.dir, "/model_", nr.states, ".txt", sep = '')
  
  # read and parse files
  model <- parseModelFile(file = file.model)
  
  # emission and overlap matrices
  emission <- model$emissions
  transition <- model$transitions
  
  ordered.marks <- c("H3K4me3", "H3K9ac", "H3K27ac", "H3K36me3", "H3K9me3", "H3K27me3")
  emission$mark.factor <- factor(x = emission$mark, levels = ordered.marks)
  
  g.emission <- ggplot(data = emission, aes(y = as.factor(state), x = mark.factor))+
    geom_tile(aes(fill = emission), col = "black")+
    theme_bw(base_size = 12)+
    scale_fill_continuous(name = "P", low = "white", high = "#4682b4", limits = c(0, 1))+
    geom_text(aes(label = round(emission, 3)), col = "black", size = 3.55)+
    xlab("Mark")+
    ylab("State")
  
  g.transition <- ggplot(data = transition, aes(y = as.factor(state.from), x = as.factor(state.to)))+
    geom_tile(aes(fill = transition), col = "black")+
    theme_bw(base_size = 12)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))+
    scale_fill_continuous(name = "P", low = "white", high = "#4682b4", limits = c(0, 1))+
    geom_text(aes(label = round(transition, 3)), col = "black", size = 3.55)+
    ylab("State (from)")+
    xlab("State (to)")
  
  return(list(g.emission = g.emission, emission = emission,
              g.transition = g.transition, transition = transition))
}




getModelLogLik <- function(dir) {
  files <- list.files(path = dir, pattern = "model_", full.names = T)
  
  model.data <- c()
  for(file in files) {
    entry <- read.csv(file = file, sep = "\t", header = F, stringsAsFactors = F)[1, ]
    states <- as.numeric(entry[1])
    log.lik <- as.numeric(entry[4])
    iterations <- as.numeric(entry[5])
    row <- c(states, iterations, log.lik)
    model.data <- rbind(model.data, row)
  }
  
  model.data <- data.frame(model.data, row.names = NULL)
  colnames(model.data) <- c("states", "iterations", "log.lik")
  return(model.data)
}




getMartGenes <- function() {
  chrs <- c("MT", "18", "17", "7", "8", "15", "6", "Y", "19", "X", "1", "16",
            "5", "11", "2", "3", "10", "9", "4", "12", "13", "14")
  
  # setup biomart database
  genes <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  attributes <- c("ensembl_gene_id", "mgi_symbol", "gene_biotype",
                  "start_position", "end_position", 
                  "chromosome_name", "strand")
  
  mart.data <- try(getBM(attributes = attributes, 
                         filters = c("chromosome_name"),
                         values = list(chrs),
                         mart = genes))
  return (mart.data)
}




getMartGenesFromGoParent <- function(go.parent.id) {
  chrs <- c("MT", "18", "17", "7", "8", "15", "6", "Y", "19", "X", "1", "16",
            "5", "11", "2", "3", "10", "9", "4", "12", "13", "14")
  
  # setup biomart database
  genes <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  attributes <- c("ensembl_gene_id", "mgi_symbol", "gene_biotype",
                  "start_position", "end_position", 
                  "chromosome_name", "strand")
  
  mart.data <- try(getBM(attributes = attributes, 
                         filters = c("chromosome_name", "go_parent_term"), 
                         values = list(chrs, go.parent.id), 
                         mart = genes))
  return (mart.data)
}




# Description:
# Given a chromosome, the procedure reads in the states returned for col2 and 
# col10, and appends them to the correct bp bin of the binarized data produced.
# Only chr1-chr19, not chrMt, chrX or chrY
getStatesChromHmm <- function(nr.states, 
                              chrom.dir) {
  chromosomes <- getChrs()
  
  states <- c()
  for(chr in chromosomes) {
    file.P <- paste(chrom.dir, "/STATEBYLINE/PC_", nr.states, 
                    "_", chr, "_statebyline.txt", sep = '')
    file.H <- paste(chrom.dir, "/STATEBYLINE/HC_", nr.states, 
                    "_", chr, "_statebyline.txt", sep = '')
    
    if(file.exists(file.P) & file.exists(file.H)) {
      states.P <- read.csv(file = file.P, stringsAsFactors = F)
      states.P <- states.P[-1, ]
      states.H <- read.csv(file = file.H, stringsAsFactors = F)
      states.H <- states.H[-1, ]
      
      states.chr <- data.frame(states.P = states.P, 
                               states.H = states.H, 
                               chr = chr, 
                               stringsAsFactors = F)
      states <- rbind(states, states.chr)
    }
  }
  
  return(states)
}




# Description:
# Given a chromosome, the procedure reads in the states returned for col2 and 
# col10, and appends them to the correct bp bin of the binarized data produced.
# Only chr1-chr19, not chrMt, chrX or chrY
getLimbStatesChromHmm <- function(nr.states, 
                                  chrom.dir) {
  chromosomes <- getChrs()
  chromosomes <- chromosomes[chromosomes != "chrMt"]
  
  states <- c()
  for(chr in chromosomes) {
    # limb
    L11.5 <- paste(chrom.dir, "/STATEBYLINE/L.11.5_", nr.states, 
                    "_", chr, "_statebyline.txt", sep = '')
    L11.5 <- read.csv(file = L11.5, stringsAsFactors = F)
    L11.5 <- L11.5[-1, ]
    
    
    L12.5 <- paste(chrom.dir, "/STATEBYLINE/L.12.5_", nr.states, 
                   "_", chr, "_statebyline.txt", sep = '')
    L12.5 <- read.csv(file = L12.5, stringsAsFactors = F)
    L12.5 <- L12.5[-1, ]
    
    
    L13.5 <- paste(chrom.dir, "/STATEBYLINE/L.13.5_", nr.states, 
                   "_", chr, "_statebyline.txt", sep = '')
    L13.5 <- read.csv(file = L13.5, stringsAsFactors = F)
    L13.5 <- L13.5[-1, ]
    
    
    L14.5 <- paste(chrom.dir, "/STATEBYLINE/L.14.5_", nr.states, 
                   "_", chr, "_statebyline.txt", sep = '')
    L14.5 <- read.csv(file = L14.5, stringsAsFactors = F)
    L14.5 <- L14.5[-1, ]
    
    
    L15.5 <- paste(chrom.dir, "/STATEBYLINE/L.15.5_", nr.states, 
                   "_", chr, "_statebyline.txt", sep = '')
    L15.5 <- read.csv(file = L15.5, stringsAsFactors = F)
    L15.5 <- L15.5[-1, ]
    
    
    # brain
    B11.5 <- paste(chrom.dir, "/STATEBYLINE/B.11.5_", nr.states, 
                   "_", chr, "_statebyline.txt", sep = '')
    B11.5 <- read.csv(file = B11.5, stringsAsFactors = F)
    B11.5 <- B11.5[-1, ]
    
    
    B12.5 <- paste(chrom.dir, "/STATEBYLINE/B.12.5_", nr.states, 
                   "_", chr, "_statebyline.txt", sep = '')
    B12.5 <- read.csv(file = B12.5, stringsAsFactors = F)
    B12.5 <- B12.5[-1, ]
    
    
    B13.5 <- paste(chrom.dir, "/STATEBYLINE/B.13.5_", nr.states, 
                   "_", chr, "_statebyline.txt", sep = '')
    B13.5 <- read.csv(file = B13.5, stringsAsFactors = F)
    B13.5 <- B13.5[-1, ]
    
    
    B14.5 <- paste(chrom.dir, "/STATEBYLINE/B.14.5_", nr.states, 
                   "_", chr, "_statebyline.txt", sep = '')
    B14.5 <- read.csv(file = B14.5, stringsAsFactors = F)
    B14.5 <- B14.5[-1, ]
    
    
    B15.5 <- paste(chrom.dir, "/STATEBYLINE/B.15.5_", nr.states, 
                   "_", chr, "_statebyline.txt", sep = '')
    B15.5 <- read.csv(file = B15.5, stringsAsFactors = F)
    B15.5 <- B15.5[-1, ]
    
    
    file.P <- paste(chrom.dir, "/STATEBYLINE/PC_", nr.states, 
                    "_", chr, "_statebyline.txt", sep = '')
    states.P <- read.csv(file = file.P, stringsAsFactors = F)
    states.P <- states.P[-1, ]
    
    
    file.H <- paste(chrom.dir, "/STATEBYLINE/HC_", nr.states, 
                    "_", chr, "_statebyline.txt", sep = '')
    states.H <- read.csv(file = file.H, stringsAsFactors = F)
    states.H <- states.H[-1, ]
    
    
    states.chr <- data.frame(states.L.11.5 = L11.5,
                             states.L.12.5 = L12.5,
                             states.L.13.5 = L13.5,
                             states.L.14.5 = L14.5,
                             states.L.15.5 = L15.5,
                             states.B.11.5 = B11.5,
                             states.B.12.5 = B12.5,
                             states.B.13.5 = B13.5,
                             states.B.14.5 = B14.5,
                             states.B.15.5 = B15.5,
                             states.P = states.P, 
                             states.H = states.H, 
                             chr = chr, stringsAsFactors = F)
    states <- rbind(states, states.chr)
  }
  
  return(states)
}




# Description:
getOverlap <- function(bin.data.genomic, 
                       bin.data.states, 
                       nr.states) {
  # merge data
  bin.data.genomic <- cbind(bin.data.genomic, bin.data.states)
  
  # labels and states to analyze
  labels <- c("cpg", "exon", "gene", "tes", "tss2kb", "tss", 
              "protein", "enhancer.id", "states.col2", "states.col10")
  nr.states <- 1:nr.states
  
  # compute % of each state in a given label
  out <- c()
  for(label in labels) {
    temp <- bin.data.genomic[bin.data.genomic[, label] == 1, ]
    col2.temp <- data.frame(table(temp$states.col2)/nrow(temp)*100, stringsAsFactors = F)
    col2.temp$label <- label
    col10.temp <- data.frame(table(temp$states.col10)/nrow(temp)*100, stringsAsFactors = F)
    col10.temp$label <- label
    
    out <- rbind(out, col2.temp)
    out <- rbind(out, col10.temp)
  }
  
  return (out)
}


