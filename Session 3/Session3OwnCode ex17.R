library(copula)
library(rugarch) #for std t density pdist(...) function
library(MASS) #for fit distrb fitdist(...) function

## Ex 17 ####

## Main Logic sub-question c and d ##

  # First fit marginal t's and store parameters (mean, shape, degrees of FD)
  # Use the fitted parameters of the univariate t's to transform data to
  # uniform (0,1) format (like an inverse CDF transform), use pdist(...) command
  # This is then the PARAMETRIC way to transform the data
  # A non parametric way would be to work with the empirical cdf like function
  
  # Calculate starting values of Omega matrix for t-copula, with Kendall's Tau.
  # We will transform Kendall's Tau values to Omega matrix values with formulas
  # of slide deck.
  
  # When determining the Omega matrix this way, we need to check for negative
  # eigenvalues. If negative eigenvalues, we replace them with a small positive
  # number, and reconstruct the correct Omega matrix via an eigenvalue-eigenvector
  # decomposition.
  
  # Finally we use the fitCopula function to find optimal Omega matrix values and
  # degrees of freedom value. The Omega matrix values calculated with Kendall's Tau
  # formula is a starting point for the (pseudo) MLE routine.

## Loading data
D <- read.csv(file.choose(),dec=",",sep=";",row.names=1)

## 17a ##
cov(D)
# AUD biggest volatility

## 17b ##
# IMPORTANT! fitdist doesn't work when the 2nd argument is a dataframe with
# a single variable, which looks very similar to a numerical vector in the R
# environment. Don't use D["JPY"] as input as this results in a dara frame
# Use D[,1] instead, which gives a numerical vector and also a valid result

EST.JPY = as.numeric(fitdist(dist="std",D[,1],"t")$par)
EST.AUD = as.numeric(fitdist(dist="std",D[,2],"t")$par)
EST.INR = as.numeric(fitdist(dist="std",D[,3],"t")$par)
# INR has lower degrees of freedom

## 17c ##

# Kendall's Tau values
KTJA = cor(D[,1],D[,2],method="kendall") #kendall's tau JPY vs AUD
KTJI = cor(D[,1],D[,3],method="kendall") #kendall's tau JPY vs INR
KTAI = cor(D[,2],D[,3],method="kendall") #kendall's tau AUD vs INR

# Omega matrix values (by transforming Kendall's tau's)
OMJA = (pi/2)*sin(KTJA) #omega matrix value JPY vs AUD
OMJI = (pi/2)*sin(KTJI) #omega matrix value JPY vs INR
OMAI = (pi/2)*sin(KTAI) #omega matrix value AUD vs INR

OM = matrix(c(1,OMJA,OMJI,OMJA,1,OMAI,OMJI,OMAI,1),ncol=3,byrow=T)

eigen(OM)
# no negative eigenvalues. So Omega matrix is OK!

## 17d ##
# transforming data in a parametric way
DU=D

# TIP: ignore the naming suggestions of the pdist function
# it goes: 3rd argument mean, 4th argument shape, 5th arugment degrees of fd
DU[,1] = pdist(dist="std",D[,1],EST.JPY[1],EST.JPY[2],EST.JPY[3])  
DU[,2] = pdist(dist="std",D[,2],EST.AUD[1],EST.AUD[2],EST.AUD[3])
DU[,3] = pdist(dist="std",D[,3],EST.INR[1],EST.INR[2],EST.INR[3])

# transforming the data in a non parametric way 
# TIP: use rank(.) function
n = nrow(D)
DUnp = matrix(c(rank(D[,1])/(n+1),rank(D[,2])/(n+1),rank(D[,3])/(n+1)),ncol=3,byrow=F)
DUnp


# fitting the copula and perfoming MLE routine
tcop = tCopula(c(OMJA,OMJI,OMAI),3,dispstr = "un",df=3)

fitMLEpar = fitCopula(tcop,DU,method="ml",start=c(OMJA,OMJI,OMAI,3))
fitMLEpar

fitMLEnonpar = fitCopula(tcop,DUnp,method="ml",start=c(OMJA,OMJI,OMAI,3))
fitMLEnonpar

ffrank = fitCopula(copula = frankCopula(3, dim = 3),data = DU, method = "ml")
ffrank



