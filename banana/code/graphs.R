library(ggplot2)
library(ggpubr)
library(latex2exp)

# colors
ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette="Set1")

# load aggregated results
means = read.csv("~/Documents/aMTM-simulations/banana/results/means.csv")
ses = read.csv("~/Documents/aMTM-simulations/banana/results/ses.csv")

# data
algos = list("None", "AM", "ASWAM", "RAM")
cols = c("None"="#E69F00", "AM"="#56B4E9", "ASWAM"="#009E73", "RAM"="#F0E442")
which = c(0, 1, 2, 3)
Ks = c(1:10)

# plot functions
plot_metric = function(metric, displayed=metric){
   plt = ggplot() + ggtitle(displayed) + labs(y="") + scale_x_continuous(
      "K", 
      breaks=c(2, 4, 6, 8, 10), 
      labels=c(2, 4, 6, 8, 10), 
      limits=c(1, 10)
   ) 
   for(i in which){
      m = means[[metric]][means$adapt==i]
      s = ses[[metric]][means$adapt==i]
      df = data.frame(K=Ks, m=m, l=m-s, u=m+s, Adaptation=algos[[i+1]])
      plt = plt + geom_ribbon(
            data=df, 
            aes(ymin=l, ymax=u, x=K, fill=Adaptation),
            alpha=0.2
         ) + geom_line(
            data=df,
            aes(x=K, y=m, colour=Adaptation)
         ) + scale_colour_manual(values=cols) + scale_fill_manual(values=cols) 
      
   }
   plt
   return(plt)
}

empty_plot = ggplot() + theme_minimal() 

no_xticks = theme(
   axis.title.x=element_blank(), 
   axis.text.x=element_blank(),
   axis.ticks.x=element_blank()
)


# create all plots
plots = list(
   plot_metric("tv", "(a) TV Distance"),
   plot_metric("ess", "(b) mESS"),
   plot_metric("msjd", "(c) MSJD"),
   plot_metric("time.elapsed", "(d) CPU"),
   plot_metric("esscpu.elapsed", "(e) mESS/CPU"),
   plot_metric("msjdcpu.elapsed", "(f) MSJD/CPU"),
   empty_plot,
   plot_metric("essnbeval", "(g) mESS/NbEval"),
   plot_metric("msjdnbeval", "(h) MSJD/NbEval")
)


# merge all plots
plts = ggpubr::ggarrange(
   plotlist=plots,
   nrow=3, ncol=3,
   common.legend=TRUE, legend="bottom",
   align="hv"
)

# save plot
ggsave(
   "~/Documents/aMTM-simulations/banana/figs/banana_results.pdf",
   plts,
   width=10, height=7, device=cairo_pdf
)


