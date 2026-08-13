#-------------------------------------------------------------

dunnett = function(k,zobs)
{
# This function concerns the distribution of k RVs obtained
# from comparing k treatments with a control.
# The standardised variates Z1, ... , Zk are assumed to
# follow a multivariate normal distribution, with each pair
# having correlation 0.5.
# The function returns the probability that the maximum
# of Z1, ... , Zk exceeds zobs.
# In the computation the Zi are treated as differences
# between sample means on each of the k treatments and
# the sample mean on the control arm.

# Set up grid points to integrate over mc = mu_hat for
# the control arm, which we take as N(0,0.5)

# Calculate Simpson's rule weights

 r=16
 out=s.rule(r,0,-50,50)
 mc=out[[2]]*sqrt(0.5)
 ww=out[[3]]*sqrt(0.5)

# Check
#
# sum(ww*dnorm(mc,0,sqrt(0.5)))

# Now compute the probability that the maximum mu_hat for
# the k treatments exceeds mc by at least zobs

 sum=0
 for(ii in 1:length(mc))
 {
  sum=sum+ww[ii]*dnorm(mc[ii],0,sqrt(0.5))*
                (1-pnorm((mc[ii]+zobs)*sqrt(2))^k)
 }
 epsilon=10^(-10)
 pdunnett=(sum <= epsilon)*epsilon +
          (sum >= 1-epsilon)*(1-epsilon) +
          (sum > epsilon & sum < 1-epsilon)*sum
 return(pdunnett)
}

#-------------------------------------------------------------

find.dunnett.c = function(k,alpha)
{
 # Finds value of k such that
 #
 #   P(Z_obs > k) = alpha
 #
 # where Z_obs is the maximum of k N(0,1) RVS with
 # correlation 0.5, as in Dunnett's test.

 clo=qnorm(1-alpha)
 flo=dunnett(k,clo)
 chi=clo*1.1
 fhi=dunnett(k,chi)
 if(fhi > alpha)
 {
  while(fhi > alpha)
  {
   clo=chi
   flo=fhi
   chi=chi*1.1
   fhi=dunnett(k,chi)
  }
 }

 while(flo-fhi > 10^(-5))
 {
  cmid=0.5*(clo+chi)
  fmid=dunnett(k,cmid)
  if(fmid > alpha)
  {
   clo=cmid
   flo=fmid
  } else
  {
   chi=cmid
   fhi=fmid
  }
 }

 croot=clo+(chi-clo)*(alpha-flo)/(fhi-flo)

 froot=dunnett(k,croot)

 return(croot)
}

#-------------------------------------------------------------