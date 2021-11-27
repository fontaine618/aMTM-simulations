# Setup
setwd("~/Documents/aMTM-simulations/banana")
source('code/problems.R')

# get sample
N = 1e3
set.seed(1)
sample = iid(N, parms)
colors = regions(sample)$regions

# plot
library(latex2exp)
library("ggplot2")
library("GGally")
data = data.frame(sample, Region=colors)
data$Region = as.factor(data$Region)
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#000000", "#0072B2", "#D55E00", "#CC79A7")
lowerFn = function(data, mapping){
    ggplot(data=data, mapping=mapping) + geom_point(alpha=0.3) +
        theme(panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank())
}
g = ggpairs(
    data, 
    columns=1:5, 
    upper=list(continuous=lowerFn, mapping=aes(colour=Region)),
    diag="blank",
    lower=list(continuous=lowerFn, mapping=aes(colour=Region)),
    legend=6,
    columnLabels=c("x[1]", "x[2]", "x[3]", "x[4]", "x[5]"),
    labeller = "label_parsed"
)
g = g + scale_colour_manual(values=cbPalette)

ggsave(
    "~/Documents/aMTM-simulations/banana/figs/banana_samples.pdf",
    g,
    width=6, height=5, device=cairo_pdf
)