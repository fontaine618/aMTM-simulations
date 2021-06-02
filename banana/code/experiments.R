library(data.table)
pdes = list("banana" = CJ(N=c(1e5)))
ades = list(
   aMTM = CJ(
      K=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20), 
      adapt=c(0, 1, 2, 3), 
      gamma=c(.7),
      accrate=c(.5),
      weight=c(0),
      proposal=c(2),
      burnin=c(1/11)
   )
)
