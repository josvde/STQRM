
library(copula)

## Ex 13a ####

# we will write a 3rd to default cumulative default function, with correlation
# as input, and then run an optim routine to find correlation where P=45%

# loading hazard rates
H <- read.csv(file.choose(),sep=";",dec=",",header=T)
H <- H$adapt2/10000

# 3rd to default cum prob calculation
P3DT <- function(p){  #correlation function
C <- normalCopula(p,125,"ex")
nsim=1000
RDC <- rCopula(nsim,C) # random draw from copula
RDDT <- apply(RDC,1,qexp,H) # transform random copula draws into default times
RDDT <- t(RDDT)
V3DT <- apply(RDDT,1,quantile,(3-1)/(125-1)) # vector with 3rd to default times
return(abs(ecdf(V3DT)(5)-0.45))}

o <- optim(0.7,P3DT,method = "L-BFGS-B",lower=0,upper=1)
o$par #answer is correlation = 0.7


## Ex 13b ####
# now we can no longer work with ecdf short cut
# we need to calculate the number of cases where
# conditions are met by hand.

nsim=1000

P13b <- function(corr){  #correlation function
  C <- normalCopula(corr,125,"ex")
  RDC <- rCopula(nsim,C) # random draw from copula
  RDDT <- apply(RDC,1,qexp,H) # transform random copula draws into default times
  RDDT <- t(RDDT)
  
  V6DT <- apply(RDDT,1,quantile,(6-1)/(125-1)) # vector with 6th to default times
  V12DT <- apply(RDDT,1,quantile,(12-1)/(125-1)) # vector with 12th to default times
  
  P = sum(V6DT<=5&V12DT>=5)/nsim
  
return(100*abs(P-0.1))}

## optim routine doesn't work!
#o <- optim(0.5,P13b,method = "L-BFGS-B",lower=0,upper=1)
#corr = o$par #answer is correlation = 0.7
#P13b(corr)

## regular optimization via pre defined grid
corrV = seq(0.60,0.75,by=0.005) # correlation grid
corrV = matrix(corrV,ncol=1)
optimR = apply(corrV,1,P13b)
optimCorr = corrV[optimR==min(optimR)]
optimCorr # optimization also yield 0.7 as optimal correlation value

rm(list=setdiff(ls(), "H"))

## Ex 14a ####

## loading hazard rates
H <- read.csv(file.choose(),sep=";",dec=",",header=T)
H <- H$adapt2/10000
## Main logic in a nutshell

  # First we need simulated default times of the entire portfolio,
  # then we will make a sub selection of the simulated default time
  # based on the tranche we are interested in.
  # We will also account for a possible partial first and possible partial last
  # credit, to make sure the tranche %'s match the credits we take into account
  # Each default corresponds to a N*(1-R) cost/loss where N is a percentage expressing 
  # the weight of 1 credit in the entire CDO/credit basket value.

## Premium calculation formula intuition:

  # We say that the value V = 0 to apply a AEPV equivalence principle.
  # The costs are the defaults needed to be covered for each defaulted credit.
  # This equals N*(1-R), because the client will receive N*R after a default,
  # and we, the seller of the CDS, cover them for the incurred losses.
  # We, the seller, are also interested in the AEPV of the entire tranche we 
  # are insuring. The tranche value (E(t)) is a sort of variable endowment. 
  # Where we discount all possible values continuously to the present.
  # We determine the one off premium (as a percentage) such that the AEPV 
  # of our losses incurred (i.e. credit losses reimbursed to client) equals the 
  # AEPV of the tranche times p. This will make sure the risk neutral value
  # of the CDS (V) is 0. In other words p is the % of the AEPV of the tranche 
  # that I need to charge to cover the AEPV of the costs (insured credit defaults).
  
## Parameter values used in the formula

  # E(t), the value of the tranche over time, will be split up in k parts E_i(t)
  # each being either 1/N, where N is total number of credits in the CDO, or 0.
  # In summary:
  # when t<Ti: Ei(t) = Ni = 1/N = 1/125
  # when t>=Ti: Ei(t) = 0 because credit defaulted
  # So the integral we need to calculated can be expressed a the sum of all Ei(t)
  # multiplied by (1-e^-rTi)/r (See notes on paper...).
  # E_i(t) behaves like a stocks price, where we continuously discount the price
  # at each point in time. So the discount factor is continuously changing.
  # In contrast with the losses, where we don't need to do this. 
  # Because they occur in 1 specific point in time, so only need to be discounted
  # that point in time.

  # Each credit loss that we need to cover for our client equals (1/125)*(1-R).
  # But for the value of the tranche itself, each credit has a (pre-loss) value of 1/125!

## Extra detail

  # p (one off premium) is expressed as a percentage of the AEPV of the tranche, needed
  # to cover the AEPV of the losses.
  # However, the AEPV of the tranche and the losses both are expressed as a % of the 
  # ENTIRE CDO, which we set to 1 for the purpose of this exercise.
  
  # In the solution code of Emanuel pay-off equals costs (as used in the explanation above)
  # and premium equals the value of the tranche (as used in explanation above).


## Simulation of sorted default times

p = 0.8 # correlation
M = 5 # maturity
r = 0.0299 # annual r
C <- normalCopula(p,125,"ex")
nsim=50000
RDC <- rCopula(nsim,C) # random draw from copula
RDDT <- apply(RDC,1,qexp,H) # transform random copula draws into default times
RDDT <- t(RDDT)
RDDTs <- apply(RDDT,1,sort) # sorted default times
RDDTs <- t(RDDTs)
D <- RDDTs # shorter variable naming...
rm(RDC,C,RDDT,RDDTs)

## Selection credit numbers according to tranche attachment and detachment points

Ni = 1/125 # entire CDO notional value = 1 & all credits have equals weight
Ni1mR = (1/125)*0.6 # cover amount we need to reimburse the client for each default
ap = 0.12 #attachment point
dp = 0.22 #detachment point

## first relevant credit calculation
# we round up because 0 after comma or even the tiniest decimal value after the comma
# indicates that the next credit in line needs to be taken into account.
# E.g. Attachment point 10% would indicates a theoretical starting credit of 12,4.
# This means we don't need the 12th credit in line, only 60% of the 13th credit, 
# and 100% of all credits above that, until we reach the detachment point.
# A theoretical first credit of 12.01, results in 99% of credit 13.
# Exactly 12 would result in 100% of credit 13.

FRC = ceiling(ap/Ni1mR+0.000001) 

## last relevant credit calculation
# This is almost the same reasoning as for the first relevant credit.
# But we don't add a small number this time because if the theoretical last credit value
# gives a whole number, we pick that credit number as the last one we need to take into account
# And we need 100% of that credit, so there will then be no fraction needed of the next credit in line.
# E.g. a theoretical last credit of 30.4 means 100% of the 30th credit in line, but only
# 40% of the 31st credit in line. 30.01 would be 100% and 1% respectively
# A value of 30 exactly, means 100% of the 30th credit and 0% of the 31st credit.

LRC = ceiling(dp/Ni1mR)

## Calculation of FRC and LRC fraction values
# For the FRC fraction we need 1 minus decimal value of theoretical first credit
# For the LRC fraction we just need the decimal value of last first credit

FRCf = FRC - ap/Ni1mR #FRC will always be higher then theoretical FRC
LRCf = dp/Ni1mR - (LRC-1) #using LRC-1 gives us the right fraction value

## Make selection of default time matrix for credits relevant to tranche
DS <- D[,FRC:LRC]

## Benefit and cost calculation according to default times

# Option 1.
# All credits mature -> cost = 0 and benefits will be maximal value
# but we don't care about benefit value as one-off premium will be 0 
# for that simulation because incurred costs are 0.

# Option 2.
# All credits default.
# Costs (Ni1mR = Ni*(1-R) = (1/125)*(1-R) ) for each credit that need to 
# be discounted according to default times. 
# Value of the tranche (when not defaulted = Ei(t) = 1/125 for each credit, 0 otherwise) 
# will need to be be discounted differently (via the integral expression just as a stock is discounted)

# Option 3. Partial defaults
# Same as option 2. but now non-defaulted credits need to be discounted for the entire lifetime of the CDS
# instead of for the simulated default times, which can be higher than the maturity.

# Option 2 & 3 can be quickly calculated together, by setting default times higher then maturity
# equal to maturity

# cost function
Cf <- function(DT){
  if(DT[1]>M){ # option 1
    res = 0
  }else{
    DTm = head(DT,length(DT)-1) #main default times without LRC (that needs a fraction)
    res=Ni1mR*sum(exp(-r*DTm[DTm<M]))
    if(tail(DT,1)<M) {
      res = res + LRCf*Ni1mR*exp(-r*tail(DT,1))}
  }
  return(res)}

# tranche value function
Bf <- function(DT){
  if(DT[1]>M){
    res = (Ni/r)*(LRC-FRC+LRCf)*(1-exp(-r*M)) # arbitrary value (as one-off premium for this simulation will be 0 because of 0 was also cost)
  }else{
    DT[DT>M] = M #set matured credit default times higher than M equal to M
    res = (Ni/r)*sum(1-exp(-r*head(DT,length(DT)-1)))
    res = res + (Ni/r)*LRCf*(1-exp(-r*tail(DT,1)))
    }
return(res)}

Cres <- apply(DS,1,Cf)
Bres <- apply(DS,1,Bf)

mean(Cres)
mean(Bres)
OneOffP = mean(Cres/Bres)
OneOffP

# mean(payoff) #in absolute terms (as a % of entire invested notional) = 20,15 bp
# [1] 0.00202006

# mean(premium)
# [1] 0.7708996

# mean(payoff/premium)# runprice = 30,7 bp (0.307%)
# [1] 0.003057376

## Ex 14b ####

# put everything from 14a in a function and optimize

OP <- function(p){
M = 5 # maturity
r = 0.0299 # annual r
C <- normalCopula(p,125,"ex")
nsim=1000
RDC <- rCopula(nsim,C) # random draw from copula
RDDT <- apply(RDC,1,qexp,H) # transform random copula draws into default times
RDDT <- t(RDDT)
RDDTs <- apply(RDDT,1,sort) # sorted default times
RDDTs <- t(RDDTs)
D <- RDDTs # shorter variable naming...
rm(RDC,C,RDDT,RDDTs)

## Selection credit numbers according to tranche attachment and detachment points

Ni = 1/125 # entire CDO notional value = 1 & all credits have equals weight
Ni1mR = (1/125)*0.6 # cover amount we need to reimburse the client for each default
ap = 0.12 #attachment point
dp = 0.22 #detachment point

FRC = ceiling(ap/Ni1mR+0.000001) 
LRC = ceiling(dp/Ni1mR)
FRCf = FRC - ap/Ni1mR #FRC will always be higher then theoretical FRC
LRCf = dp/Ni1mR - (LRC-1) #using LRC-1 gives us the right fraction value

## Make selection of default time matrix for credits relevant to tranche
DS <- D[,FRC:LRC]

Cres <- apply(DS,1,Cf)
Bres <- apply(DS,1,Bf)

mean(Cres)
mean(Bres)
OneOffP = mean(Cres/Bres)

return(abs(OneOffP-0.06))}

Ores <- optim(0.6,OP,method="L-BFGS-B",lower=0,upper=1)

# optimal correl +/- 0.57

















