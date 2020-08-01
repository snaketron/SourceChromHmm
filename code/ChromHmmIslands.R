# requirements
source("code/util/GeneralUtilities.R")
source("code/util/ChromHmmUtilities.R")
require("doParallel")
require("foreach")
require("doMC")
registerDoMC(cores = 12)



# Description:
# Given the states data, this procedure combines contiguous states into islands.
getIslands <- function(keys, chr) {
  temp.keys <- keys
  island.id <- numeric(length = length(keys))
  i <- 1
  count <- 1
  while(i < length(keys)) {
    j <- keys == keys[i]
    k <- min(which(j == F))-1
    if(is.infinite(k)) {
      k <- length(keys)
    }
    
    island.id[i:k] <- count
    keys[i:k] <- NA
    count <- count + 1
    i <- k + 1
    cat(i, "\n")
  }
  
  return (data.frame(chr = chr, 
                     island.id = island.id,
                     keys = temp.keys,
                     stringsAsFactors = F))
}



#### Mutations islands HC+PC ####
# states data
state.data <- get(load(file = "output/states.15.RData"))
state.data$states.M <- paste(state.data$states.P, state.data$states.H, sep = '.') # mutations

# find mutation islands
chrs <- unique(state.data$chr)
islands <- (foreach(i = 1:length(chrs)) %dopar% getIslands(keys = state.data$states.M[state.data$chr == chrs[i]], chr = chrs[i]))
islands <- do.call(rbind, islands)
state.data[, "island.M"] <- paste(islands$chr, islands$keys, islands$island.id, sep = '.')

state.data <- state.data[, c("chr", "bp.start", "bp.stop", "states.M", "island.M")]
colnames(state.data) <- c("chr", "bp.start", "bp.stop", "state", "island")

states.min <- aggregate(bp.start~chr+state+island, data = state.data, FUN = min)
states.max <- aggregate(bp.stop~chr+state+island, data = state.data, FUN = max)
islands <- merge(x = states.min, y = states.max, by = c("chr", "state", "island"))

islands <- islands[order(islands$bp.start, decreasing = F), ]
islands <- islands[order(islands$chr, decreasing = F), ]
islands$bp.mid <- islands$bp.start+(islands$bp.stop-islands$bp.start)/2
x <- do.call(rbind, strsplit(x = islands$state, split = "\\."))
islands$state.P <- x[, 1]
islands$state.H <- x[, 2]
save(islands, file = "output/chromhmm_islands/mutation.islands.15.RData")
rm(x, state.data, chrs, islands, states.min, states.max)



#### Mutations islands HC+PC+Limb+Brain ####
# states data
state.data <- get(load(file = "output/states.LB.15.RData"))
state.data$states.M <- paste(state.data$states.P, state.data$states.H, sep = '.') # mutations

# find mutation islands
chrs <- unique(state.data$chr)
islands <- (foreach(i = 1:length(chrs)) %dopar% getIslands(keys = state.data$states.M[state.data$chr == chrs[i]], chr = chrs[i]))
islands <- do.call(rbind, islands)
state.data[, "island.M"] <- paste(islands$chr, islands$keys, islands$island.id, sep = '.')

state.data <- state.data[, c("chr", "bp.start", "bp.stop", "states.M", "island.M")]
colnames(state.data) <- c("chr", "bp.start", "bp.stop", "state", "island")

states.min <- aggregate(bp.start~chr+state+island, data = state.data, FUN = min)
states.max <- aggregate(bp.stop~chr+state+island, data = state.data, FUN = max)
islands <- merge(x = states.min, y = states.max, by = c("chr", "state", "island"))

islands <- islands[order(islands$bp.start, decreasing = F), ]
islands <- islands[order(islands$chr, decreasing = F), ]
islands$bp.mid <- islands$bp.start+(islands$bp.stop-islands$bp.start)/2
x <- do.call(rbind, strsplit(x = islands$state, split = "\\."))
islands$state.P <- x[, 1]
islands$state.H <- x[, 2]
save(islands, file = "output/chromhmm_islands/mutation.islands.LB.15.RData")
rm(x, state.data, chrs, islands, states.min, states.max)






#### State islands HC+PC ####
# states data
state.data <- get(load(file = "output/states.15.RData"))

chrs <- unique(state.data$chr)
p <- (foreach(i = 1:length(chrs)) %dopar% getIslands(keys = state.data$states.P[state.data$chr == chrs[i]], chr = chrs[i]))
h <- (foreach(i = 1:length(chrs)) %dopar% getIslands(keys = state.data$states.H[state.data$chr == chrs[i]], chr = chrs[i]))
p <- do.call(rbind, p)
h <- do.call(rbind, h)

state.data[, "island.PC"] <- paste(p$chr, p$keys, p$island.id, sep = '.')
state.data[, "island.HC"] <- paste(h$chr, h$keys, h$island.id, sep = '.')

p <- state.data[, c("chr", "bp.start", "bp.stop", "states.P", "island.PC")]
colnames(p) <- c("chr", "bp.start", "bp.stop", "state", "island")
h <- state.data[, c("chr", "bp.start", "bp.stop", "states.H", "island.HC")]
colnames(h) <- c("chr", "bp.start", "bp.stop", "state", "island")

states.min <- aggregate(bp.start~chr+state+island, data = p, FUN = min)
states.max <- aggregate(bp.stop~chr+state+island, data = p, FUN = max)
islands.PC <- merge(x = states.min, y = states.max, by = c("chr", "state", "island"))

states.min <- aggregate(bp.start~chr+state+island, data = h, FUN = min)
states.max <- aggregate(bp.stop~chr+state+island, data = h, FUN = max)
islands.HC <- merge(x = states.min, y = states.max, by = c("chr", "state", "island"))

save(islands.PC, file = "output/chromhmm_islands/PC.islands.15.RData")
save(islands.HC, file = "output/chromhmm_islands/HC.islands.15.RData")
rm(p, h, states.max, states.min, islands.PC, islands.HC, chrs, state.data)


#### State islands HC+PC+Limb+Brain ####
state.data <- get(load(file = "output/states.LB.15.RData"))

chrs <- unique(state.data$chr)
p <- (foreach(i = 1:length(chrs)) %dopar% getIslands(keys = state.data$states.P[state.data$chr == chrs[i]], chr = chrs[i]))
h <- (foreach(i = 1:length(chrs)) %dopar% getIslands(keys = state.data$states.H[state.data$chr == chrs[i]], chr = chrs[i]))
p <- do.call(rbind, p)
h <- do.call(rbind, h)

state.data[, "island.PC"] <- paste(p$chr, p$keys, p$island.id, sep = '.')
state.data[, "island.HC"] <- paste(h$chr, h$keys, h$island.id, sep = '.')

p <- state.data[, c("chr", "bp.start", "bp.stop", "states.P", "island.PC")]
colnames(p) <- c("chr", "bp.start", "bp.stop", "state", "island")
h <- state.data[, c("chr", "bp.start", "bp.stop", "states.H", "island.HC")]
colnames(h) <- c("chr", "bp.start", "bp.stop", "state", "island")

states.min <- aggregate(bp.start~chr+state+island, data = p, FUN = min)
states.max <- aggregate(bp.stop~chr+state+island, data = p, FUN = max)
islands.PC <- merge(x = states.min, y = states.max, by = c("chr", "state", "island"))

states.min <- aggregate(bp.start~chr+state+island, data = h, FUN = min)
states.max <- aggregate(bp.stop~chr+state+island, data = h, FUN = max)
islands.HC <- merge(x = states.min, y = states.max, by = c("chr", "state", "island"))

save(islands.PC, file = "output/chromhmm_islands/PC.islands.LB.15.RData")
save(islands.HC, file = "output/chromhmm_islands/HC.islands.LB.15.RData")
rm(p, h, states.max, states.min, islands.PC, islands.HC, chrs, state.data)



##### states data #####
# do for each combination 15 vs. LB.15

# state.data <- get(load(file = "output/states.15.RData"))
state.data <- get(load(file = "output/states.LB.15.RData"))

chrs <- unique(state.data$chr)
keys <- colnames(state.data)[4:15]
for(key in keys) {
  islands <- (foreach(i = 1:length(chrs)) %dopar% getIslands(
    keys = state.data[state.data$chr == chrs[i], key], chr = chrs[i]))
  islands <- do.call(rbind, islands)
  state.data[, paste("island", key, sep = '.')] <- paste(state.data$chr, islands$island.id, sep = '.')
  cat(key, "\n")
}


keys <- colnames(state.data)[which(regexpr(pattern = "island\\.", text = colnames(state.data)) != -1)]
for(key in keys) {
  state.data[, key] <- paste(state.data[, "chr"], state.data[, key], sep = '.')
  cat(key, "\n")
}
# save(state.data, file = "output/chromhmm_islands/states.15.islands.RData")
save(state.data, file = "output/chromhmm_islands/states.LB.15.islands.RData")
rm(islands, key, keys)


