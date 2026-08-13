# Example 2: Survival data
#
# Two-arm trial, new treatment vs control
#
# Observations follow a proportional hazards model with
# log hazard ratio theta, where
#
#    hazard ratio > 1 and theta > 0
#
# imply the new treatment is superior
#
# To test H_0: theta =< 0 vs theta > 0 with one-sided type I error
# probability alpha=0.025 and power 1-beta=0.8 at
#
#   theta = delta = 0.5 = log(1.65)

# Load functions to be called later

source("~/My files/deming_dec_2025_demos/standards.R")
source("~/My files/deming_dec_2025_demos/ersptest.R")
source("~/My files/deming_dec_2025_demos/findroot.R")

alpha=0.025
beta=0.2
delta=0.5

# Design a fixed sample size trial

ifix=(qnorm(1-alpha)+qnorm(1-beta))^2/delta^2
ifix

# ifix=31.40

# Design a group sequential trial with 5 analyses
#
# Error spending design, non-binding futility boundary, rho-family
# error spending function, rho=2 for efficacy and futility boundaries
#
# Call function find_rfac_ersp1b to find the inflation factor

r=16         #  Gives high accuracy for numerical computations
na=5         #  Number of analyses
alpha=0.025  #  Type I error rate
beta=0.2     #  Type II error rate
rho=2        #  Determines rate of error spending

rfac=find_rfac_ersp1b(r,na,alpha,beta,rho)
rfac

# Required inflation factor is rfac=1.133 (the same to 3dp as for beta=0.1)

imax=rfac*ifix
imax

# imax=35.58

# Inspect this design with cumulative information levels as planned

iobs=imax*c(1:5)/5
iobs

# Define cumulative error probabilities cume
#
# cume[2;1:na]  cumulative type 1 error probabilities under theta=0
# cume[1;1:na]  cumulative type 2 error probabilities under theta=delta
#
# If we observe data with iobs[na] > imax, we cap type I error spent at
# alpha and the type II error spent at beta
#
# If we observe data with iobs[na] < imax, we set the final type I error
# to be alpha, so cume[2,na=alpha] 

rho=2
cume=array(0,c(2,na))
cume[2,1:na]=pmin(alpha*(iobs/imax)^rho,alpha)
cume[2,na]=alpha
cume[1,1:na]=pmin(beta*(iobs/imax)^rho,beta)
cume

out=ersp1b(r,na,iobs,delta,cume)
out  # Boundary, type I error probability, type II error probability

zbdy=out[[1]]
zbdy

plot(c(1:na),zbdy[1,],col="blue",ylim=c(-2,4),xlim=c(0.5,5.2),pch=20,
     main="Error spending boundary",ylab="Z",xlab="Analysis")
points(c(1:na),zbdy[2,],col="red",pch=20)

# Expected information on termination

theta.vec=seq(-0.5,1.5,0.04)
ntheta=length(theta.vec)
einf.vec=rep(NA,ntheta)
for(ii in c(1:ntheta))
{
  theta=theta.vec[ii]
  einf.vec[ii]=gst1(r,na,iobs,zbdy,theta)[[4]]
}
plot(theta.vec,einf.vec,ylim=c(0,40),pch=20,col="blue",
     main="E(Inf on termination)",ylab="E(Inf)",xlab="theta")
lines(theta.vec,rep(31.4,ntheta),col="red")

# Running the trial with observed information

iobs=c(5.43,12.58,21.11,30.55,33.28)
rho=2
cume=array(0,c(2,na))
cume[2,1:na]=pmin(alpha*(iobs/imax)^rho,alpha)
cume[2,na]=alpha
cume[1,1:na]=pmin(beta*(iobs/imax)^rho,beta)
cume

out=ersp1b(r,na,iobs,delta,cume)
zbdy=out[[1]]
zbdy

# With a binding futility boundary

rfac=find_rfac_ersp1c(r,na,alpha,beta,rho)
rfac

imax=rfac*ifix
imax

cume=array(0,c(2,na))
cume[2,1:na]=pmin(alpha*(iobs/imax)^rho,alpha)
cume[2,na]=alpha
cume[1,1:na]=pmin(beta*(iobs/imax)^rho,beta)
cume

out=ersp1c(r,na,iobs,delta,cume)
zbdy=out[[1]]
zbdy

# Back to a non-binding futility boundary

rfac=find_rfac_ersp1b(r,na,alpha,beta,rho)
rfac

imax=rfac*ifix
imax

# Analysis of Cox model output

iobs=c(4.11,10.89,19.23,28.10,30.96)
rho=2
cume=array(0,c(2,na))
cume[2,1:na]=pmin(alpha*(iobs/imax)^rho,alpha)
cume[2,na]=alpha
cume[1,1:na]=pmin(beta*(iobs/imax)^rho,beta)
cume

out=ersp1b(r,na,iobs,delta,cume)
zbdy=out[[1]]
zbdy
