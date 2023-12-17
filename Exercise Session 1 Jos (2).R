### exercise session 1

## 1a ####

# loading data

AAPLdata = read.csv(file.choose(),sep=",",dec=".")
VZdata = read.csv(file.choose(),sep=",",dec=".")

AAPLdata = AAPLdata$Adj.Close
VZdata = VZdata$Adj.Close

# compute DV portfolio based in final price (last entry in data vectors!)
AAPLcv = tail(AAPLdata,1) # closing value for dollar value portfolio
VZcv = tail(VZdata,1) # closing value for dollar value portfolio

DV = 1000*VZcv+500*AAPLcv
w1 = 1000*VZcv/DV # Verizon weight
w2 = 500*AAPLcv/DV # Apple weight
W = c(w1,w2)

rm(AAPLcv,VZcv,w1,w2)

# transform prices in DV into % returns
end = length(AAPLdata)
AAPLr = (AAPLdata[2:end]-AAPLdata[1:(end-1)])/AAPLdata[1:(end-1)]
AAPLr = AAPLr*100
VZr = (VZdata[2:end]-VZdata[1:(end-1)])/VZdata[1:(end-1)]
VZr = VZr*100

# portfolio return
mu = c(mean(VZr),mean(AAPLr))
mup = t(W)%*%mu

# portfolio variance
data = matrix(c(VZr,AAPLr),ncol=2)
COV = cov(data)
COV
pv = t(W)%*%COV%*%W

# 1 day VaR
VaR1a = -DV*qnorm(0.01,mup,sqrt(pv))/100
VaR1a

# 1 day ES
n=10000
a=0.01
VaR1av= rep(0,n)

for(i in 1:n){
  VaR1av[i]=-DV*qnorm(a*(i/n),mup,sqrt(pv))/100
}
ES1a = mean(VaR1av)
ES1a

rm(a,i,n,end)

## 1b ####

# USE QUANTILE FUNCTION ON UNSORTED DATA TO DETERMINE HISTORICAL VAR

portData = W[1]*data[,1]+W[2]*data[,2] # vector with portfolio returns
VaR1b = quantile(portData,0.01)*(-DV)/100 
VaR1b

portDataDV = portData*DV/100 # transform % portfolio returns into changes in portfolio DV
ES1b = portDataDV[portDataDV <= - VaR1b] # select all observed DV changes below or equal to VAR
ES1b = -mean(ES1b) # take average to get ES
ES1b

rm(COv,mup,pv,VaR1a,VaR1b,VaR1av,VaR1bmu,qdp,ES1a,ES1b,mu,portData,portDataDV,COV)

## 2a ####

library('MASS')
library('mnormt')

# transform data into non % returns, to make sure scaling matrix is in correct units
# to calculate VAR

dataRNP = data/100 #data in non percentage returns
rm(data)

# fit bivariate T
fit = cov.trob(dataRNP,nu=3,cor=TRUE)
muv3 = fit$center
COVv3 = 3*fit$cov

# results
muv3
COVv3

## 2b ####

# General info (on 2b & 2c)

  # There are 2 way to do the calibration: 1. On the returns and 2. On the dollar values.
  
  # I will calibrate on returns (not expressed in %!), to make sure no mistake are made
  # w.r.t. to the scaling matrix units.
  
  # Because you can't calculate quantiles of a generalized T (with a scaling lambda and mu),
  # one VERY important step to remember is: calculate a quantile for a classical 
  # univariate T (via qt command). And then transform this quantile into a quantile of 
  # a generalized univariate T, by using the LAMBDA OF THE PORTFOLIO!
  # Not THE STD DEV OF THE PORTFOLIO!

# Finding optimal nu (degrees of freedom)
nuV = seq(3.20,3.23,by=0.001) #nu vector grid
LLV = rep(0,length(nuV)) # log likelihood vector
for(i in 1:length(nuV)){
  fit = cov.trob(dataRNP,nu=nuV[i])
  LLV[i] = sum(log(dmt(dataRNP,fit$center,fit$cov,df=nuV[i])))
}
MaxLLV = max(LLV)
nuOpt = nuV[LLV==MaxLLV]
nuOpt # 3.21

# resulting fitted mu and lambda
fitOpt = cov.trob(dataRNP,nu=nuOpt)
muOpt = fitOpt$center
SOpt = fitOpt$cov

muOpt
SOpt

## 2c ####

# calculate portfolio parameters
muPort = muOpt%*%W
lambdaPort = t(W)%*%SOpt%*%W
lambdaPort = sqrt(lambdaPort) # DO NOT FORGET TO TAKE SQUARE ROOT!
muPort
lambdaPort

# VAR via quantile transformation
Q99 = qt(0.01,df = nuOpt)
VaR2cPer = Q99*lambdaPort+muPort
VaR2cPer
VaR2cDV = -DV*VaR2cPer
VaR2cDV #3580 OK
