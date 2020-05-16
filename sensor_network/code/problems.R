#############################################
# Generating the problem 2, 22, 6, 25, 26 38!

n_sensors = 6
n_known = 2
n_unknown = 4
known_sensors = c(F,F,F,F,T,T)
set.seed(6)
X = round(matrix(
   c(
      0.57,0.91,
      0.10,0.37,
      0.24,0.14,
      0.85,0.04,
      0.50,0.30,
      0.30,0.70
      ),
   nrow=n_sensors,
   ncol=2, T
),2)
y = as.matrix(dist(X, diag=F, upper=F))

probs = exp(-y^2 / (2*0.5^2))
w = apply(probs, 1:2, function(p) rbinom(1,1,p))
w[upper.tri(w)] = 0

w = diag(n_sensors)
w[c(4,5,6), 1] = 1
w[c(3,4), 2] = 1
w[c(5,6), 3] = 1


y = y*w
y[y==0]=NA
y[upper.tri(y)] = NA

segs = data.frame(i=NA, j=NA, x0=NA, y0=NA, x1=NA, y1=NA)


for(i in seq(2,n_sensors)){
   for(j in seq(1, i-1)){
      if(w[i,j] == 1){
         segs = rbind(segs, c(i,j,X[i,1],X[i,2],X[j,1],X[j,2]))
      } 
   }
}


# Define the target
Xk = as.vector(t(X[known_sensors,]))
Xu = as.vector(t(X[!known_sensors,]))

parms=list(Xk=Xk, y=y, w=w, sigy=0.02, sigw=0.3, sigx=9, 
           ns=n_sensors, nu=n_unknown, nk=n_known)

logp = function(x, parms){
   K = length(x)/(2*parms$nu)
   xu= x
   x = cbind(x, matrix(parms$Xk, K, parms$nk*2, byrow=T))
   X = array(x, dim=c(K, 2, parms$ns))
   
   out = sapply(seq(K), function(k){
      d = as.matrix(dist(t(X[k,,]), diag=T, upper=T))
      llky = -sum((parms$y - d)^2 / (2*parms$sigy^2), na.rm=T)
      llkw = -sum(parms$w*d^2 / (4*parms$sigw^2), na.rm=T)
      llkw1 = sum((1-parms$w) * log(1-exp(-d^2 / (2*parms$sigw^2))), na.rm=T)
      prior = -exp(sum((xu[k, ]-0.5)^2) / (2*parms$sigx))
      #prior = -sum(abs(xu[k, ]) > 1) * 1e+10
      return(llky + llkw + llkw1 + prior)
   })
   return(out)
}


# precompile for faster evaluation
target <- compiler::cmpfun(logp)

# prepare data and function
data = list(
   parms = parms,
   target = logp
)
fun = function(data, job, ...){
   list()
}


