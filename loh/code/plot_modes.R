library(ggplot2)
library(egg)
library(latex2exp)
library(directlabels)
library(RColorBrewer)
library(ggpubr)
library(rmutil)
# Setup
setwd("~/Documents/aMTM-simulations/loh")
df = read.table("data/BarrettsLOH.dat", header=F)
colnames(df) = c("x", "n")
df$p = df$x/df$n
# where to evaluate the mixture
m = median(df$n)
# modes
m1 = list(eta=0.903, pi1=0.228, pi2=0.708, gamma=3.54)
m2 = list(eta=0.078, pi1=0.832, pi2=0.230, gamma=-18.51)
# evaluate the mixtures
x = seq(0, m)
bin1 = dbinom(x, m, m1$pi1)
bbin1 = dbetabinom(x, m, m1$pi2, 2*(1+exp(m1$gamma))/exp(m1$gamma))
mix1 = m1$eta * bin1 + (1-m1$eta) * bbin1
d1 = data.frame(bin=bin1, bbin=bbin1, mix=mix1, x=x/m)
bin2 = dbinom(x, m, m2$pi1)
bbin2 = dbetabinom(x, m, m2$pi2, 2*(1+exp(m2$gamma))/exp(m2$gamma))
mix2 = m2$eta * bin2 + (1-m2$eta) * bbin2
d2 = data.frame(bin=bin2, bbin=bbin2, mix=mix2, x=x/m)

plot_curves = function(d, tit){
    plt = ggplot(data=df, aes(p))
    plt = plt + geom_histogram(aes(y = stat(count) / sum(count)), breaks=seq(0, 16)/16)
    plt = plt + labs(y = "Frequency", x = "Proportion") + ggtitle(tit)
    plt = plt + geom_line(data=d, aes(x=x, y=bin, linetype="Binomial"))
    plt = plt + geom_line(data=d, aes(x=x, y=bbin, linetype="BetaBinomial"))
    plt = plt + geom_line(data=d, aes(x=x, y=mix, linetype="Mixture"))
    plt = plt + scale_linetype_manual(values=c(
        Binomial=2, BetaBinomial=3, Mixture=1
        ), guide=guide_legend(fill = NULL, colour = NULL))
    plt = plt + theme(
        legend.title = element_blank(),
        legend.position = c(.50, .95),
        legend.justification = c("top"),
        legend.margin = margin(0, 6, 6, 6)
    )
    plt
}

plt1 = plot_curves(d1, "Mode 1 (0.903, 0.228, 0.708, 3.54)")
plt2 = plot_curves(d2, "Mode 2 (0.078, 0.832, 0.230, -18.51)")

plt = ggarrange(
    plt1, plt2,
    nrow=1
)

ggsave(
    "~/Documents/aMTM-simulations/loh/figs/hist_modes.pdf",
    plt,
    width=10, height=3, device=cairo_pdf
)

