#############################
## Formats ChromHMM output ##
#############################

# source
source("code/util/ChromHmmUtilities.R")
source("code/util/GeneralUtilities.R")

# const
nr.states <- 15




#  PC vs HC
chrom.dir <- "output/chromhmm"

# now construct the states
header <- get(load(file = "output/header.RData"))

state.data <- getStatesChromHmm(nr.states = nr.states,
                                chrom.dir = chrom.dir)
state.data$chr <- NULL
header <- cbind(header, state.data)
state.data <- header
save(state.data, file = paste0("output/states.", nr.states, ".RData"))





# PC vs HC vs LIMB vs BRAIN
chrom.dir <- "output/chromhmm_limb_brain/"

# now construct the states
header <- get(load(file = "output/header.RData"))

state.data <- getLimbStatesChromHmm(nr.states = nr.states,
                                    chrom.dir = chrom.dir)
state.data$chr <- NULL
header <- cbind(header, state.data)
state.data <- header
save(state.data, file = paste0("output/states.LB.", nr.states, ".RData"))
