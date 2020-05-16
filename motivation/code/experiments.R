library(data.table)
pdes = list("motivation" = CJ(N=c(1e4)))
ades = list(
   IID = CJ(),
   Metropolis = CJ(scales=list(list(400)), adapt=list(0), sorted=FALSE),
   AM = CJ(scales=list(list(400)), adapt=list(1), sorted=FALSE),
   MTM = CJ(scales=list(list(400, 100, 9)), adapt=list(0), sorted=FALSE),
   aMTM = CJ(scales=list(list(400, 36, 36)), adapt=list(3), sorted=FALSE)
)
