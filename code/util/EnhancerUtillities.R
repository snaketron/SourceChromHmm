# Description:
# For each row of the subject, it tests if it has any overlap 
# (parial, relaxed), with any row of the query.
getPartialOverlapCount <- function(query, 
                                   subject, 
                                   key, 
                                   with.relaxed = T) {
  out <- c()
  chrs <- unique(query$chr)
  for(chr in chrs) {
    temp.subject <- subject[subject$chr == chr, ]
    temp.query <- query[query$chr == chr, ]
    temp.query[, paste(key, "count", sep = '.')] <- 0
    if(with.relaxed) {
      temp.query[, paste("relaxed", key, "count", sep = '.')] <- 0
    }
    
    q <- IRanges(start = temp.query$bp.start, end = temp.query$bp.stop)
    s <- IRanges(start = temp.subject$bp.start, end = temp.subject$bp.stop)
    
    # partial overlap
    m <- findOverlaps(query = q, subject = s, type = "any")
    if(length(m@from) != 0) {
      m.count <- data.frame(table(m@from), stringsAsFactors = F)
      colnames(m.count) <- c("query.id", "count")
      m.count$query.id <- as.numeric(as.character(m.count$query.id))
      temp.query[m.count$query.id, paste(key, "count", sep = '.')] <- m.count$count
    }
    
    
    # relaxed overlap
    if(with.relaxed) {
      mr <- findOverlaps(query = q, subject = s, type = "any", maxgap = 200)
      if(length(mr@from) != 0) {
        mr.count <- data.frame(table(mr@from), stringsAsFactors = F)
        colnames(mr.count) <- c("query.id", "count")
        mr.count$query.id <- as.numeric(as.character(mr.count$query.id))
        temp.query[mr.count$query.id, paste("relaxed", key, "count", sep = '.')] <- mr.count$count
      }
    }
    
    # collect
    out <- rbind(out, temp.query)
  }
  
  return (out)
}



# Description:
# For each row of the subject, it tests if it has any overlap 
# (parial, relaxed), with any row of the query.
getLabelledOverlap <- function(subject, 
                               query, 
                               key, 
                               with.relaxed) {
  
  out <- c()
  chrs <- unique(subject$chr)
  for(chr in chrs) {
    temp.subject <- subject[subject$chr == chr, ]
    temp.subject[, key] <- NA
    
    if(with.relaxed) {
      temp.subject[, paste("relaxed", key, sep = '.')] <- NA
    }
    
    temp.query <- query[query$chr == chr, ]
    s <- IRanges(start = temp.subject$bp.start, end = temp.subject$bp.stop)
    q <- IRanges(start = temp.query$bp.start, end = temp.query$bp.stop)
    
    
    # partial overlap
    m <- findOverlaps(query = q, subject = s, type = "any")
    m <- data.frame(as.matrix(m))
    if(nrow(m) != 0) {
      m$ID <- temp.query$ID[m$queryHits]
      m <- aggregate(ID~subjectHits, data = m, FUN = paste, collapse = ';')
      temp.subject[unique(m$subjectHits), key] <- m$ID
    }
    
    if(with.relaxed) {
      # relaxed overlap
      relax.m <- findOverlaps(query = q, subject = s, type = "any", maxgap = 200)
      relax.m <- data.frame(as.matrix(relax.m))
      if(nrow(relax.m) != 0) {
        relax.m$ID <- temp.query$ID[relax.m$queryHits]
        relax.m <- aggregate(ID~subjectHits, data = m, FUN = paste, collapse = ';')
        temp.subject[unique(relax.m$subjectHits), paste("relaxed", key, sep = '.')] <- relax.m$ID
      }
    }
    
    # collect
    out <- rbind(out, temp.subject)
  }
  
  return (out)
}





# Description:
# For each row of the subject, it tests if it has any overlap 
# (parial, relaxed), with any row of the query.
getPartialOverlapEnhancerCount <- function(query, 
                                           subject, 
                                           key) {
  out <- c()
  chrs <- unique(query$chr)
  for(chr in chrs) {
    temp.subject <- subject[subject$chr == chr, ]
    temp.query <- query[query$chr == chr, ]
    temp.query[, paste(key, "sum", sep = '.')] <- 0
    
    q <- IRanges(start = temp.query$bp.start, end = temp.query$bp.stop)
    s <- IRanges(start = temp.subject$bp.start, end = temp.subject$bp.stop)
    
    # partial overlap
    m <- findOverlaps(query = q, subject = s, type = "any")
    if(length(m@from) != 0) {
      m <- data.frame(from = m@from, to = m@to, to.count = temp.subject$Count[m@to])
      m <- aggregate(to.count~from, data = m, FUN = sum)
      temp.query[m$from, paste(key, "sum", sep = '.')] <- m$to.count
    }
    
    # collect
    out <- rbind(out, temp.query)
  }
  
  return (out)
}





# Description:
# For each row of the subject, it tests if it has any overlap 
# (parial, relaxed), with any row of the query.
getPartialOverlapLabel <- function(query, 
                                   subject, 
                                   key) {
  out <- c()
  chrs <- unique(query$chr)
  for(chr in chrs) {
    temp.subject <- subject[subject$chr == chr, ]
    temp.query <- query[query$chr == chr, ]
    temp.query[, paste(key, "overlap", sep = '.')] <- NA
    
    q <- IRanges(start = temp.query$bp.start, end = temp.query$bp.stop)
    s <- IRanges(start = temp.subject$bp.start, end = temp.subject$bp.stop)
    
    # partial overlap
    m <- findOverlaps(query = q, subject = s, type = "any")
    if(length(m@to) != 0) {
      m <- data.frame(from = m@from, to = m@to)
      m <- m[duplicated(m) == F, ]
      m$to <- temp.subject$ID[m$to]
      m <- aggregate(to~from, data = m, FUN = base::paste, collapse = ';')
      temp.query[m$from, paste(key, "overlap", sep = '.')] <- m$to
    }
    
    # collect
    out <- rbind(out, temp.query)
  }
  
  return (out)
}




# Description:
getOrderedTadGenes <- function(bp.mid, 
                               tad.id, 
                               tad.data, 
                               genes.data) {
  genes <- ''
  distances <- ''
  if(!is.na(tad.id)) {
    genes <- tad.data$genes[tad.data$ID == tad.id]
    genes <- unlist(strsplit(genes, split = ';'))
    
    if(length(genes) != 0) {
      distances <- bp.mid - genes.data$tss[genes.data$Gene.name %in% genes]
      
      # sort
      genes <- genes[order(abs(distances))]
      distances <- distances[order(abs(distances))]
      genes <- paste(genes, collapse = ";")
      distances <- paste(distances, collapse = ";")
    }
  }
  
  return(list(genes = genes, distances = distances))
}



# Description:
# For each row of the subject, it tests if it has any overlap 
# (parial, relaxed), with any row of the query.
getLabelledOverlap <- function(subject, 
                               query, 
                               key, 
                               with.relaxed) {
  
  out <- c()
  chrs <- unique(subject$chr)
  for(chr in chrs) {
    temp.subject <- subject[subject$chr == chr, ]
    temp.subject[, key] <- NA
    
    if(with.relaxed) {
      temp.subject[, paste("relaxed", key, sep = '.')] <- NA
    }
    
    temp.query <- query[query$chr == chr, ]
    s <- IRanges(start = temp.subject$bp.start, end = temp.subject$bp.stop)
    q <- IRanges(start = temp.query$bp.start, end = temp.query$bp.stop)
    
    
    # partial overlap
    m <- findOverlaps(query = q, subject = s, type = "any")
    m <- data.frame(as.matrix(m))
    if(nrow(m) != 0) {
      m$ID <- temp.query$ID[m$queryHits]
      m <- aggregate(ID~subjectHits, data = m, FUN = paste, collapse = ';')
      temp.subject[unique(m$subjectHits), key] <- m$ID
    }
    
    if(with.relaxed) {
      # relaxed overlap
      relax.m <- findOverlaps(query = q, subject = s, type = "any", maxgap = 200)
      relax.m <- data.frame(as.matrix(relax.m))
      if(nrow(relax.m) != 0) {
        relax.m$ID <- temp.query$ID[relax.m$queryHits]
        relax.m <- aggregate(ID~subjectHits, data = m, FUN = paste, collapse = ';')
        temp.subject[unique(relax.m$subjectHits), paste("relaxed", key, sep = '.')] <- relax.m$ID
      }
    }
    
    # collect
    out <- rbind(out, temp.subject)
  }
  
  return (out)
}



