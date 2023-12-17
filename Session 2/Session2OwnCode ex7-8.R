## Set up 
library('MASS')
library('mnormt')

## 7a ####
data = read.csv(file.choose(),sep=";",dec=",")
data = data[,2:6]

# Profile Likelihood optimization
nuV = seq(3.55,3.6,by=0.001) #nu vector grid
LLV = rep(0,length(nuV)) # log likelihood vector
for(i in 1:length(nuV)){
  fit = cov.trob(data,nu=nuV[i])
  LLV[i] = sum(log(dmt(data,fit$center,fit$cov,df=nuV[i])))
}
MaxLLV = max(LLV)
nuOpt = nuV[LLV==MaxLLV]
nuOpt = round(nuOpt,digits=2)
nuOpt #3.57

# resulting fitted mu and lambda
fitOpt = cov.trob(data,nu=nuOpt)
muOpt = fitOpt$center
SOpt = fitOpt$cov

muOpt
SOpt

## 7b ####

# compute weights (DV of positions and DV of total portfolio are given)
W = c(0.2,0.1,0.2,0.2,0.3)
DV=5000000
# portfolio t parameters
muPort = muOpt%*%W
lambdaPort = t(W)%*%SOpt%*%W
lambdaPort = sqrt(lambdaPort) # DO NOT FORGET TO TAKE SQUARE ROOT!
muPort
lambdaPort

# ES calculation with VAR function and VAR discrete averaging
fVaR <- function(a){
Q99 = qt(a,df = nuOpt)
VaR2cPer = Q99*lambdaPort+muPort
VaR2cPer
VaR2cDV = -DV*VaR2cPer
return(VaR2cDV)}

aGrid = seq(0.0000001,0.02,by=0.0000001)
VaRv = rep(0,length(aGrid)) #var vector
for(i in 1:length(aGrid)){
  VaRv[i]= fVaR(aGrid[i])
}
ES98 = mean(VaRv)
ES98 #381561 OK

rm(fit,aGrid,ES98,DV,LLV,MaxLLV,nuV,VaRv,W,fVaR,i,fitOpt,lambdaPort,muPort,SOpt,muOpt,nuOpt)

## 8a ####

# General framework 8a:

  # Important for this sub-question is transformation from the scale matrix to 
  # the cholesky decomposition matrix A and also the transformation back to the scale matrix 
  # Formulas -> see LAST slide multivariate statistics slide deck!
  # Cholesky decomposition for any matrix in R via "chol()" command
  
  # Step 1 in this exercise is to Write a log likelihood function that can be used
  # in optim command. This function takes in 6 parameters:
  # 3 for cholesky's A, 2 for mu vector and 1 for nu
  # Also see last slide multivariate statistics slide deck for formulas for these values,
  # where d = 2
  # In this function Cholesky's A is constructed from the input vector P.
  # This cholesky A matrix is then used to recconstruct the scale matrix S needed for
  # the LL optimization

  # Step 2 is the optimization itself.
  # Starting values for A are given by the cholesky decomposition on sample covariance matrix
  # For mu, we use the sample mean return.
  # TIP: if optim doesn't converge --> set starting values equal to the results of a cov.trob 
  # routine. Meaning, determine A from fitted scale matrix, and use fitted mu and nu.
  # optim command specs: USE method = BGFS and set Hessian = TRUE
  # When optim still doesn't converge after using appropriate start vales, boundaries can be
  # added by setting method to L-BGFS-B

  # Step 3 (optional): compare cov.trob fitted parameters and max LL value with those 
  # of the optim fit routine 

# Data is sub set of data set ex. 7
data = data[,1:2]


# Step 1 Log Likelihood function (with parameter vector P with dimension 6)
LL <- function(P){
  A = matrix(P[1],P[2],0,P[3],ncol=2,byrow=T) # reconstruct cholesky
  S = t(A)%*%A # construct scale matrix
  mu = c(P[4],P[5])
  nu = P[6]
  LLsum = sum(log(dmt(data,mu,S,nu)))
  return(LLsum)}






