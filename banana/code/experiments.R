library(data.table)
pdes = list("banana" = CJ(N=c(1e5)))
ades = list(
   aMTM = CJ(
      K=1:10, 
      adapt=c(0, 1, 2, 3), 
      gamma=c(.7),
      accrate=c(.2),
      weight=c(0),
      proposal=c(2),
      burnin=c(1/2),
      local=F,
      global=F
   )
)
