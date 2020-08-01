# State transition matrices 
# 1) state changes
# 2) island changes


# mutation states
state.data <- get(load(file = "output/states.15.RData"))
state.data$states.P <- factor(x = state.data$states.P, levels = 1:15)
state.data$states.H <- factor(x = state.data$states.H, levels = 1:15)
m <- as.data.frame.matrix(table(state.data$states.P, state.data$states.H))
save(m, file = "output/mutation_matrix/MutationStates.15.RData")
rm(m, state.data)



# mutation islands
islands <- get(load(file = "output/chromhmm_islands/mutation.islands.15.RData"))
islands$states.P <- factor(x = islands$state.P, levels = 1:15)
islands$states.H <- factor(x = islands$state.H, levels = 1:15)
m <- as.data.frame.matrix(table(islands$states.P, islands$states.H))
save(m, file = "output/mutation_matrix/MutationIslands.15.RData")
rm(m, islands)

