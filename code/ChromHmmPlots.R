####################
## ChromHMM plots ##
####################


# required libs
source("code/util/ChromHmmUtilities.R")
source("code/util/GeneralUtilities.R")
require("ggplot2")
require("gridExtra")
require("reshape2")
# require("extrafont")
# loadfonts(device = "pdf")

# global variables
nr.states <- 15

# chrom.dir <- "output/chromhmm"
# key <- ''
# states <- get(load(file = "output/states.15.RData"))


chrom.dir <- "output/chromhmm_limb_brain/"
key <- 'LB'
states <- get(load(file = "output/states.LB.15.RData"))



# Plot Likelihood #
model.data <- getModelLogLik(dir = chrom.dir)
model.data <- model.data[order(model.data$states, decreasing = F), ]
model.data <- model.data[model.data$states > 3, ]
model.data$AIC <- -2*model.data$log.lik + 2*(model.data$states * (model.data$states-1) + model.data$states*7)


# Log likelihood
g1 <- ggplot(data = model.data)+
  geom_point(aes(x = states, y = log.lik))+
  theme_bw(base_size = 12)+
  scale_x_continuous(name = "Number of states", 
                     breaks = model.data$states, 
                     labels = model.data$states)+
  ylab(label = "Log likelihood")


# Akaike
g2 <- ggplot(data = model.data)+
  geom_point(aes(x = states, y = AIC))+
  theme_bw(base_size = 12)+
  scale_x_continuous(name = "Number of states", 
                     breaks = model.data$states, 
                     labels = model.data$states)+
  ylab(label = "Akaike information criterion (AIC)")



# load data
et <- getEmissionsAndTransitions(nr.states = nr.states, chrom.dir = chrom.dir)
rm(state.data)
names.p <- names(table(states$states.P))
p <- as.numeric(table(states$states.P))
h <- as.numeric(table(states$states.H))
counts.df <- data.frame(state = as.numeric(names.p), count.p = p, count.h = h,
                        prop.p = p/nrow(states)*100, prop.h = h/nrow(states)*100,
                        count.ratio = p/h,
                        count.ratio.label = paste(round(x = p/h, digits = 3),
                                                  " (", round(p, digits = 2), 
                                                  "/", round(h, digits = 2), ")", 
                                                  sep = ''),
                        cscr = paste(' ', paste(rep(x = ' ', times = 11), collapse = ''), 
                                     ' ', collapse = ''))


g3 <- ggplot(data = counts.df, aes(y = as.factor(state), x = "Ratio (PC/HC)"))+
  geom_tile(fill = "white", col = "black")+
  theme_bw(base_size = 12)+
  geom_text(aes(label = count.ratio.label), col = "black", size = 3.55)+
  xlab("")+
  ylab("State")




# save tiff
ggsave(filename = paste0("output/chromhmm_plots/LogLik_", key, nr.states, ".tiff"),
       plot = g1, device = "tiff", width = 5, height = 4, dpi = 300)


# save tiff
ggsave(filename = paste0("output/chromhmm_plots/AIC_", key, nr.states, ".tiff"),
       plot = g2, device = "tiff", width = 5, height = 4, dpi = 300)


# save tiff
tiff(file = paste0("output/chromhmm_plots/E_", key, nr.states, ".tiff"),
     width = 9, height = 7, units = "in", res = 300)
grid.arrange(et$g.emission, g3, layout_matrix = matrix(data = c(1, 1, 2), 
                                                       nrow = 1, byrow = T))
dev.off()

# save tiff
tiff(file = paste0("output/chromhmm_plots/T_", key, nr.states, ".tiff"), 
    width = 7, height = 6, units = "in", res = 300)
et$g.transition
dev.off()
