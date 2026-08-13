#-------------------------------------------------------------------------------

mesh = function (r) 
{
 # This function creates a mesh of points for use in integration by Simpson's
 # rule. The mesh is appropriate for an integrand which is equal to a N(0,1)
 # density multiplied by a smooth function. 
 # In applying Simpson's rule, each pair of successive mesh points defines an
 # interval. The integrand is evaluated at the ends of this interval and at its
 # midpoint and these values are weighted in the ratio 1:2:1.
 #
 # Input:  r     Used to determine the number of mesh points.
 #               If the integrand is precisely the N(0,1) density, a value
 #               of r = 16 should evaluate probabilities within 10**-6 of
 #               their true values. Approximately one decimal place of accuracy
 #               is lost in reducing r by a factor of two.
 #               Similar accuracy is to be expected if the integrand is a N(0,1)
 #               density multipled by a smooth function taking values in the
 #               range 0 to 1. 
 #
 # Output: mesh(1:6r-2)    Mesh of points for use with Simpson's rule.

 aa = -3-4*log(r/(1:(r-1)))
 bb= -3+3*(0:(4*r))/(2*r)
 cc = 3+4*log(r/((r-1):1))
 mesh = c(aa, bb, cc)

 mesh
}

#-------------------------------------------------------------------------------

gst1 = function(r,na,inf,zbdy,theta)
{
 # This function calculates properties of a group sequential boundary.
 #
 # Input:  r              Used to determine the number of mesh points (see
 #                        the function mesh)
 #         na             Number of analyses
 #         inf(1:na)      Information at the sequence of analyses
 #         zbdy(1:2,1:na) Boundary values on the Z-scale
 #                          zbdy(1,1:na) = lower boundary
 #                          zbdy(2,1:na) = upper boundary
 #         theta          Treatment effect
 #
 #                        The Z-statistic at analysis k, Z_k, is distributed
 #                        as N(theta sqrt(inf(k)),1), the joint distribution
 #                        of the Z_ks is multivariate normal, and increments
 #                        in S_k = Z_k * sqrt(inf(k)) are independent.
 # 
 #
 # Output: pstop(1:3,1:na)  pstop(1,k)=P(Trial stops by crossing the lower
 #                                       boundary at analysis k)
 #                          pstop(2,k)=P(Trial stops by crossing the upper
 #                                       boundary at analysis k)
 #                          pstop(3,k)=pstop(1,k)+pstop(2,k)
 #         pu               P(Exit by upper boundary)
 #         pl               P(Exit by lower boundary)
 #         einf             E(information on termination)

 pl=0
 pu=0
 pstop=array(0,c(3,na))
 einf=0

 k1=0
 inf1=0
 m1=1
 z1=array(0, c(1))
 h1=array(1, c(1))
 k2=1

 while(k2 <= na)
 {
  inf2=inf[k2]
  dinf=inf2-inf1
  z2var=1-inf1/inf2
  z2sd=sqrt(z2var)
  ez2=z1*sqrt(inf1/inf2)+theta*(dinf/sqrt(inf2))
  pupper=sum((1-pnorm((zbdy[2,k2]-ez2)/z2sd))*h1)
  pu=pu+pupper
  plower=sum(pnorm((zbdy[1,k2]-ez2)/z2sd)*h1)
  pl=pl+plower
  pstop[2,k2]=pstop[2,k2]+pupper
  pstop[1,k2]=pstop[1,k2]+plower
  einf=einf+dinf*sum(h1)
  if(k2 < na)
  {
   zmesh=theta*sqrt(inf2)+mesh(r)                             # Mesh, analysis k2
   zmeshb=zmesh[zmesh >= zbdy[1,k2] & zmesh <= zbdy[2,k2]]    # Modify mesh
   if(min(zmesh) < zbdy[1,k2]) {zmeshb=c(zbdy[1,k2],zmeshb)}
   if(max(zmesh) > zbdy[2,k2]) {zmeshb=c(zmeshb,zbdy[2,k2])}
   len=length(zmeshb)
   m2=2*len-1
   z2=array(0,c(m2))
   for(i in 1:len) {z2[2*i-1]=zmeshb[i]}
   if(len > 1)
   {
    for(i in 1:(len-1)) {z2[2*i]=0.5*(z2[2*i-1]+z2[2*i+1])}
   }
   m2=2*len-1                                            # Simpson's rule weights
   w2=array(0, c(m2))
   w2[1]=0
   if(m2>1)
   {
    w2[1]=(z2[3]-z2[1])/6
    w2[m2]=(z2[m2]-z2[m2-2])/6
    if(m2>3)
    {
     for(i2 in seq(3,m2-2,2)) {w2[i2]=(z2[i2+2]-z2[i2-2])/6}
    }
    for(i2 in seq(2,m2-1,2)) {w2[i2]=(z2[i2+1]-z2[i2-1])*4/6}
   }
   h2=array(0, c(m2))
   for(i2 in 1:m2)
   {
    h2[i2]=sum(dnorm(z2[i2]-ez2,sd=z2sd)*w2[i2]*h1)
   }
   k1=k2
   inf1=inf2
   m1=m2
   z1=z2
   h1=h2
  }
  k2=k2+1
 }
 pstop[3,]=pstop[1,]+pstop[2,]

 list(
      pstop,
      pu,
      pl,
      einf
      )
}

#-------------------------------------------------------------------------------