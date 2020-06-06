library(ggplot2)
library(egg)
library(latex2exp)
library(directlabels)
library(RColorBrewer)
library(ggpubr)
# Setup
setwd("~/Documents/aMTM-simulations/loh")
source('code/problems.R')
N = 1e4
d = 2
res = list()
target = logp


p = 4
K = 7
theta = matrix(
    c(runif(K), runif(K), runif(7), runif(K, -30, 30)),
    K, p, F 
)

logp(theta, parms)


eta_gam = data.frame(
    eta=c(0.9, 0.6, 0.25, 0.1),
    gam=c(5, 0, 0 ,-15)
)



# colors
lvls=c(-5, 5) - 100
ggplot <- function(...) ggplot2::ggplot(...) + scale_fill_gradient(low = "#000000", high = "#000000")
# formatting wrapper
frm = function(x, i) format(round(x, i), nsmall=i)
# return contour plot
plot_contour = function(e_g, res=100, levels=lvls){
    eta = e_g[1]
    gam = e_g[2]
    pis = seq(0.001, 0.999, length.out=res)
    theta = as.matrix(expand.grid(eta = eta, pi1 = pis, pi2 = pis, gamma = gam))
    vals = logp(theta, parms) - 100
    vals_mat = matrix(vals, res, res)
    df = data.frame(theta)
    df$logp = vals
    
    tit = paste("\u03b7 = ", frm(eta, 2), ", \u03b3 = ", frm(gam, 0), sep="")
    
    plt = ggplot()
    plt = plt + geom_contour(aes(x=pi1, y=pi2, z=logp), data=df, breaks=levels, 
                             color="black")
    plt = plt + labs(y = TeX("$\\pi_2$"), x = TeX("$\\pi_1$"))
    
    plt = plt + annotate("label", x=1, y=1, label=tit, hjust=1, vjust=1, 
                         fill=alpha(c("#EEEEEE"),0.8), label.size=0)
    plt = plt + xlim(0, 1) + ylim(0, 1)
    plt = plt + stat_contour(
        data=df, 
        aes(x=pi1, y=pi2, z=logp, colour=..level..), 
        breaks=levels) + 
        scale_color_gradient(low = "#000000", high = "#000000")
    plt = direct.label(plt, "bottom.pieces")
    plt
}



plts = apply(eta_gam, 1, plot_contour)

plt1 = plts[[1]] + theme()
plt2 = plts[[2]] + theme(
    axis.title.y=element_blank(), 
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
)
plt3 = plts[[3]] + theme(
    axis.title.y=element_blank(), 
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
)
plt4 = plts[[4]] + theme(
    axis.title.y=element_blank(), 
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
)


plt = ggarrange(
    plt1, plt2, plt3, plt4,
    nrow=1
)

plt = annotate_figure(
    plt,
    top=text_grob("LOH target (log) density contours", hjust=1.8)
    )

ggsave(
    "~/Documents/aMTM-simulations/loh/figs/contours.pdf",
    plt,
    width=10, height=3, device=cairo_pdf
)





#############################################

theta0 = c(0.078, 0.83, 0.23, -18)
theta0 = c(0.897, 0.229, 0.714, 15.661)

theta0 = c(runif(3), rnorm(1,0,10))

K = 3

log_seq = function(start, end, n){
    l = seq(log(start), log(end), length.out = n)
    exp(l)
}

sig0 = sapply(log_seq(1, 0.001, K), function(sig2){
    diag(c(0.08, 0.04, 0.06, 300)) * sig2
}, simplify="array")
#sig0[, , 1] = sig1


N = 1e4

mcmc = aMTM::aMTM(
    target = data$target,
    x0 = theta0,
    parms = data$parms,
    sig0=sig0,
    K=K,
    N=N,
    global=T,
    local=F,
    scale=F,
    accrate=0.3,
    proposal=2,
    gamma=0.5,
    adapt=3,
    burnin=0.1,
    weight=-1
)

X = mcmc$X









ids = seq(1,N,by=N/1000)
aMTM::plot.aMTM(mcmc, vars = 1:4, type='b')
#aMTM::plot.aMTM(mcmc, vars = 1:4, pairs=F)

mcmc$sel.prop
mcmc$acc.rate

mcmc$Sig
theta0

ids = seq(1,N,by=100)
plot(as.matrix(mcmc$X[ids,2:3]), col=1, xlim=0:1, ylim=0:1)