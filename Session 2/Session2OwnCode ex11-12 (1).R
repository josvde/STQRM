## Ex 11 ####
c10 = normalCopula(0.2,10,dispstr = "ex")
nsim = 100000
RDC = rCopula(nsim,c10)
u5 = pexp(5,0.01)
RDCmin = apply(RDC,1,min)
pFTD5Y = 1 - sum(RDCmin>=u5)/nsim #cum default probability first to default before 5Y.
pFTD5Y #answer = 33.2%

# Note you can also use f=ecdf(X), which assigns an empirical cdf function to f
# When you then use f(a) the cum pr of X < a is given.

## Ex 12 ####
#load hazard rates
H <- read.csv(file.choose(),sep=";",dec=",",header=F,
              col.names = c("credit number","hazard rate"))
H <- H$hazard.rate 

c12 <- normalCopula(0.3,10,"ex")
nsim = 50000
RDC <- rCopula(nsim,c12) #random draw copula
RDDT <- apply(RDC,1,qexp,H) #random draw transformed to default times
# where each column has its own hazard rate
RDDT <- t(RDDT)

V3DT <- apply(RDDT,1,quantile,2/9)# 3rd to default vector
ecdf(V3DT)(5) #answer 3.45%
