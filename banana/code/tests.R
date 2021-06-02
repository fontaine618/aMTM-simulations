logpi = function(x, p) apply(x,1,function(x) {
    -0.5 * (
        x[1]^2/(2*p$a1^2) +
        (p$B1*x[1]^2+x[2])^2 +
        x[3]^2 +
        x[4]^2/(2*p$a2^2) +
        (p$B2*x[4]^2+x[5])^2
    )
})

iid = function(n, p){
    x = matrix(rnorm(n*5), n, 5)
    x[, 1] = x[, 1] * p$a1
    x[, 2] = x[, 2] + p$B1 * (x[, 1]^2)
    x[, 4] = x[, 4] * p$a2
    x[, 5] = x[, 5] + p$B2 * (x[, 4]^2)
    return(x)
}


p = list(a1=1, a2=1, B1=3, B2=1)


sample = iid(1000000, p)

ranges = apply(sample, 2, range)

x = sample
reg1 = (x[, 4]>=0) * (x[, 5]>=5)
reg2 = (x[, 4]<0) * (x[, 5]>=5)
reg3 = (x[, 1]>=0) * (x[, 2]>=10)
reg4 = (x[, 1]<0) * (x[, 2]>=10)
reg5 = (1-reg1)*(1-reg2)*(1-reg3)*(1-reg4)
c(mean(reg1), mean(reg2), mean(reg3), mean(reg4), mean(reg5))

colors = rep(1, nrow(sample))
colors[reg1==1] = 2
colors[reg2==1] = 3
colors[reg3==1] = 4
colors[reg4==1] = 5

pairs(sample[1:10000, ], col=colors[1:10000])


