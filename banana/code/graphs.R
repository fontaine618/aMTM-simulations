library(ggplot2)
library(egg)
library(latex2exp)
library(mixtools)
# import res
load("~/Documents/aMTM-simulations/motivation/results/samples.Rdata")

# colors
ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette="Set1")

# return plot
frm = function(x, i) format(round(x, i), nsmall=i)
plot_sample = function(name, tit){
   ress = res[[name]]
   df = data.frame(ress["X"])
   df$sel = as.factor(ress["sel"]$sel)
   txt = paste("Acc. rate = ", frm(ress["results"]$results[10], 3),
               "\nMSJD = ", frm(ress["results"]$results[7], 3),
               "\nP(x\u2081>5) = ", frm(ress["results"]$results[1], 3),
               sep="")
   plt = ggplot(df, aes(x=X.1, y=X.2)) +
      geom_point(shape=19, alpha=0.2, show.legend = FALSE)
   plt = plt + labs(y = TeX("$x_2$"), x = TeX("$x_1$")) + ggtitle(tit)
   if(!is.null(ress$sig)){
      K = dim(ress$sig)[3]
      for(k in seq(K)){
         ell = ellipse(c(10, 5), ress$sig[, , k] * ress$lam[k], alpha=0.5, draw=FALSE)
         ell = data.frame(ell[1:249, ], ell[2:250, ], sel=as.factor(k))
         plt = plt + 
            geom_segment(data=ell, 
                         aes(x=X1, y=X2, xend=X1.1, yend=X2.1), 
                         show.legend = FALSE)
      }
   }
   plt = plt + annotate("label", x=30, y=20, label=txt, hjust=1, vjust=1, 
                        fill=alpha(c("#EEEEEE"),0.8), label.size=0)
   plt = plt + xlim(-5, 30) + ylim(-5, 20) 
   return(plt)
}

iid = plot_sample("IID", "(a) I.I.D.") + theme(
   axis.title.x=element_blank(), 
   axis.text.x=element_blank(),
   axis.ticks.x=element_blank()
)

metropolis = plot_sample("Metropolis", "(b) Metropolis") + theme(
   axis.title.x=element_blank(), 
   axis.text.x=element_blank(),
   axis.ticks.x=element_blank(),
   axis.title.y=element_blank(), 
   axis.text.y=element_blank(),
   axis.ticks.y=element_blank()
)

am1 = plot_sample("AM (one mode)", "(c) AM, one mode") + theme(
   axis.title.x=element_blank(), 
   axis.text.x=element_blank(),
   axis.ticks.x=element_blank(),
   axis.title.y=element_blank(), 
   axis.text.y=element_blank(),
   axis.ticks.y=element_blank()
)

am2 = plot_sample("AM (two modes)", "(d) AM, two modes")

mtm = plot_sample("MTM", "(e) MTM") + theme(
   axis.title.y=element_blank(), 
   axis.text.y=element_blank(),
   axis.ticks.y=element_blank()
)

amtm = plot_sample("aMTM", "(f) aMTM") + theme(
   axis.title.y=element_blank(), 
   axis.text.y=element_blank(),
   axis.ticks.y=element_blank()
)



plt = ggarrange(
   iid,
   metropolis,
   am1,
   am2,
   mtm, 
   amtm,
   nrow=2
)

ggsave(
   "~/Documents/aMTM-simulations/motivation/figs/samples.pdf",
   plt,
   width=10, height=5.5, device=cairo_pdf
)

   