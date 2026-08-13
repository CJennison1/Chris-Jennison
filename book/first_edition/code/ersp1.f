c=======================================================================
c
      subroutine ersp1(na,inf,delta,cume,alpha,beta,ind,r,na2,zbdy,ierr)
c
c-----------------------------------------------------------------------
c
c   ERSP1 computes a boundary for a group sequential one-sided test
c   according to the error spending approach. 
c
c   One-sided tests have some complications not seen in two-sided tests.
c   These issues and the way we deal with them are explained below.
c
c                                         C. Jennison, 24 September 1999
c
c-----------------------------------------------------------------------
c
c   INPUT
c
c      Number of analyses = na                       INTEGER
c      Information levels = inf(1),...,inf(na)       REAL*8(100)
c      Positive value of theta at which power is set
c                         = delta                    REAL*8(100)
c
c   Cumulative error probabilities are specified for exiting
c
c      upper boundary under theta=0:      cume(1;1,...,na)
c      lower boundary under theta=delta:  cume(2;1,...,na)
c
c                                                    REAL*8(2,100)
c
c   Overall total Type I and II error probabilities are also set
c
c      total Type I error prob  = alpha              REAL*8(100)
c      total Type II error prob = beta               REAL*8(100)
c
c   The test will be allowed to stop if ever this can be done with
c
c      Type I error prob = alpha   and   Type II error prob =< beta.
c
c   Indicator variable
c
c      ind = 1  for an interim analysis              INTEGER
c            2  for the final possible analysis
c
c   At the final analysis, the test will be forced to terminate by
c   convergence of the upper and lower boundaries.
c
c   Standardized statistics Z(1),...,Z(na) are assumed to follow the
c   canonical joint distribution
c
c      (Z(1),...,Z(na)) ~ multivariate normal,
c      Z(k) ~ N(theta sq(inf(k)),1),   k=1,...,na
c      cov(Z(k1),Z(k2)) = sq(inf(k1)/inf(k2)),   for k1 =< k2.
c
c   The parameter 
c
c      r                                            INTEGER
c
c   determines the number of grid points in numerical integrations, the
c   maximum number of points used being 12r-3.
c   We recommend the value r=16 to obtain probabilties to an accuracy
c   of about 10**-6.  About 1 decimal place is lost in reducing r by
c   a factor of two.  Lower accuracy may occur if increments in inf(k)
c   are small, in particular, if
c
c      (inf(k)-inf(k-1))/inf(k) < 0.01.
c
c-----------------------------------------------------------------------
c
c   NOTES
c
c   The most direct way to construct an error spending boundary is to
c   compute boundary values at analyses 1, 2, ... in turn ensuring that
c   the specified cumulative error probabilities under theta=0 and delta
c   are met at each stage.
c
c   At the final analysis input should include cume(1;na) = alpha and
c   cume(2;na) = beta.  In general, the upper and lower boundaries may
c   not meet up exactly, in which case a total Type I error of alpha is
c   attained by using the value for the upper boundary (calculating
c   probabilities under theta=0) in both upper and lower boundaries.
c
c   Several complications can arise in this approach.
c
c     1. Upper and lower boundaries may cross:
c
c     We should then terminate the study at the current analysis, using
c     the value for the upper boundary as the single critical value.
c     This gives priority to setting the total Type I error, letting
c     the power take up the effects of varying information sequences.
c
c     However, if this should happen at an interim, rather than final,
c     analysis, the total Type I error will only be cume(1,na), not the
c     desired alpha.
c
c     2. Boundaries may be at + or - infinity:
c
c     If Pr(continue to analysis na) under theta=0 is less than the
c     increment cume(1,na)-cume(1,na-1), it is not possible to spend
c     all of this Type I error at analysis na even if the test stops at
c     this analysis with rejection of H0 in all cases.  The best one
c     can do is set both boundaries at minus infinity in order to spend
c     as much additional type I error as possible.
c
c     An analogous situation occurs if Pr(continue to analysis na) under
c     theta=delta is less than cume(2,na)-cume(2,na-1) and a boundary
c     at plus infinity is then called for.
c
c   In order to resolve these issues, we have followed the approach
c   described at the end of Section 7.3.1 in the book Group Sequential
c   Methods with Applications to Clinical Trials by Jennison & Turnbull:
c
c     At each new analysis, k, we first check whether the test can be
c     brought to a conclusion with total error rates no greater than
c     alpha and beta.  If so, this is done, setting the boundary to give
c     Type I error probability exactly equal to alpha. If not, upper and
c     lower boundaries are found in the standard way to give cumulative
c     Type I and II error probabilities of cume(1;k) and cume(2;k).
c
c-----------------------------------------------------------------------
c
c   PERMISSIBLE  VALUES
c
c      0 < na =< 100
c      inf(k) > 0 for each k and inf(k) strictly increasing with k
c      delta > 0
c
c      0 < cume(1,k) =< alpha for each k
c      0 < cume(2,k) =< beta  for each k
c      both cume(1,k) and cume(2,k) strictly increasing with k
c
c      0 < alpha < 1 and  0 < beta < 1
c
c      ind = 1 or 2
c
c      r in the range 3 to 16.
c
c-----------------------------------------------------------------------
c
c   OUTPUT
c
c   Number of analyses actually used
c
c      na2
c
c   The stopping boundary is specified as
c
c      upper:  zbdy(2,1),...,zbdy(2,na)             REAL*8(100)
c      lower:  zbdy(1,1),...,zbdy(1,na)             REAL*8(100)
c
c      ierr       error indicator                   INTEGER
c
c                    =1  if an error is detected
c                    =0  if no error occurs
c
c-----------------------------------------------------------------------
c
c   SUBROUTINE and FUNCTION calls
c
c      NORCDF
c
c-----------------------------------------------------------------------
c

      integer na,r,ierr,k,k1,k2,j,m1,m2,i1,i2,lim1,lim2,nmesh,nmesh2,
     +       j1,j2,i,ii
      real*8 inf(100),cume(2,100),zbdy(2,100),theta(2),
     +       z1(2,200),z2(2,200),h1(2,200),h2(2,200),inf1,inf2,ez2(200),
     +       z2var,pu,pl,b1,b2,b3,f1,f2,f3,fstar,zmean,zmesh(100),
     +       w2(200),zrg,epslon,
     +       norcdf,f,x,v

      f(x,v)=dexp(-x*x/(2.0d0*v))/dsqrt(6.2831853d0*v)

c  Set convergence measure

      epslon=0.000001

c   Check input variables.

      ierr=0

      if(na.lt.1 .or. na.gt.100) ierr=1
      if(ierr.eq.1) goto 980

      do 1 k=1,na
1        if(inf(k).le.0) ierr=1
      if(ierr.eq.1) goto 982

      if(na.ge.2) then
         do 2 k=2,na
2           if(inf(k).le.inf(k-1)) ierr=1
         if(ierr.eq.1) goto 984
      endif

      if(delta.le.0) ierr=1
      if(ierr.eq.1) goto 986

      do 3 k=1,na
         if(cume(2,k).le.0 .or. cume(2,k).gt.beta) ierr=1
         if(cume(1,k).le.0 .or. cume(1,k).gt.alpha) ierr=1
3     continue
      if(ierr.eq.1) goto 988

      if(na.ge.2) then
         lim2=na-1
         do 4 k=1,lim2
            if(cume(2,k).le.cume(2,k-1)) ierr=1
            if(cume(1,k).le.cume(1,k-1)) ierr=1
4        continue
      endif
      if(ierr.eq.1) goto 988

      if(alpha.le.0 .or. alpha.ge.1 .or. beta.le.0 .or. beta.ge.1)
     +   ierr=1
      if(ierr.eq.1) goto 990

      if(ind.lt.1 .or. ind.gt.2) ierr=1
      if(ierr.eq.1) goto 992

      if(r.lt.1 .or. r.gt.16) ierr=1
      if(ierr.eq.1) goto 994

c  Initialise boundaries, exit probabilities, etc

      do 5 i=1,na
      do 5 j=1,2
5        zbdy(j,i)=0.0

      pl=0.0
      pu=0.0

      zrg=3.0+4.0*dlog(r*1.0d0)

c   Iterative calculations

c   At each iteration, we first find the new boundary points that meet
c   the error spending specifications.  Then, variables are updated
c   for use in numerical calculations at the next stage

c   We use separate numerical integration grids for calculations under
c   theta=0 and theta=delta.

c   For theta=0:
c
c   an integral at analysis k of e(z_k) over the sub-distribution of z_k,
c   for paths continuing to analysis k, is approximated by the sum
c
c      sum(i=1,m_k(1))  h_k(1,i_k) e(z_k(1,i_k))
c
c   For theta=delta:
c
c   an integral at analysis k of e(z_k) over the sub-distribution of z_k,
c   for paths continuing to analysis k, is approximated by the sum
c
c      sum(i=1,m_k(2))  h_k(2,i_k) e(z_k(2,i_k))
c
c   The general step is to compute the two grids of z values and
c   associated weights
c
c      z_k(1,i) and h_k(1,i),   i=1,m_k(1)  for use with theta=0
c
c   and
c
c      z_k(2,i) and h_k(2,i),   i=1,m_k(2)  for use with theta=delta.
c
c   We take the generic step to be from analysis k1 to k2=k1+1.
c
c   Notation: we write
c
c      h1     for   h_k1
c      m1     for   m_k1
c      inf1   for   inf(k1)
c      z1     for   z(k1)
c
c      h2     for   h_k2
c      m2     for   m_k2
c      inf2   for   inf(k2)
c      z2     for   z(k2)
c
c   The conditional distribution of Z2 given Z1=z1 is normal with
c
c      mean = z1*sq(inf1/inf2) + theta*((inf2-inf1)/sq(inf2))  and
c
c      variance = 1 - inf1/inf2.

c   Initialisations

      k1=0
      inf1=0.0
      m1=1
      z1(1)=0.0
      h1(1)=1.0

c  Loop step from k1 to k2=k1+1

50    continue

      k2=k1+1
      inf2=inf(k2)

      z2var=1.0-inf1/inf2

      do 60 i1=1,m1
         ez2(i1)=z1(i1)*dsqrt(inf1/inf2)
60    continue

c  Find lower boundary point, zbdy(1,k2)

c  Let
c
c     f(b) = Pr{Exit lower boundary at analysis k2}
c
c  when zbdy(1,k2)=b.
c
c  Then we must search for the value b such that
c
c     f(b) = fstar = cume(1,k2)-pl
c
c  where pl is the cumulative probability of exiting the lower boundary
c  up to analysis k1.

      fstar=cume(1,k2)-pl

c  Find a lower bound, b1, for b

c  Try the lowest z value used in integrating over z.  If this is too
c  high, set the boundary here anyway.

      b1=-zrg
      f1=0.0
      do 100 i1=1,m1
100      f1=f1+h1(i1)*norcdf((b1-ez2(i1))/dsqrt(z2var))
      if(f1.ge.fstar) then
         zbdy(1,k2)=-zrg
         goto 180
      endif

c  Find an upper bound, b2, for b

c  First try b2=0.
c  If this is too low, try the highest value used in integrating over z.
c  If this value is too low, set the boundary here anyway.

      b2=0.0
      f2=0.0
      do 110 i1=1,m1
110      f2=f2+h1(i1)*norcdf((b2-ez2(i1))/dsqrt(z2var))
      if(f2.le.fstar) then
         b2=zrg
         f2=0.0
         do 120 i1=1,m1
120         f2=f2+h1(i1)*norcdf((b2-ez2(i1))/dsqrt(z2var))
         if(f2.lt.fstar) then
            zbdy(1,k2)=zrg
            goto 180
         endif
      endif

c  Iteration to narrow the range (b1,b2).  Initially, we have
c
c     f1 < fstar < f2
c
c  and each iteration preserves this condition.

150   continue

c     print *,'150:b1,b2,f1,f2 ',b1,b2,f1,f2

      if(f2-f1 .le. epslon) goto 170

      b3=0.5*(b1+b2)
      f3=0.0
      do 160 i1=1,m1
160      f3=f3+h1(i1)*norcdf((b3-ez2(i1))/dsqrt(z2var))

      if(f3.lt.fstar) then
         b1=b3
         f1=f3
      else
         b2=b3
         f2=f3
      endif

      goto 150

170   zbdy(1,k2)=b1+(b2-b1)*(fstar-f1)/(f2-f1)

180   continue

c  Find upper boundary point, zbdy(2,k2)

c  Now, we let
c
c     f(b) = Pr{Exit upper boundary at analysis k2}
c
c  when zbdy(2,k2)=b.
c
c  Then we must search for the value b such that
c
c     f(b) = fstar = cume(2,k2)-pu
c
c  where pu is the cumulative probability of exiting the upper boundary
c  up to analysis k1.

      fstar=cume(2,k2)-pu

c  Find an upper bound, b1, for b

      b1=zrg
      f1=0.0
      do 200 i1=1,m1
200      f1=f1+h1(i1)*(1-norcdf((b1-ez2(i1))/dsqrt(z2var)))
      if(f1.ge.fstar) then
         zbdy(2,k2)=zrg
         goto 280
      endif

c  Find a lower bound, b2, for b

      b2=0.0
      f2=0.0
      do 210 i1=1,m1
210      f2=f2+h1(i1)*(1-norcdf((b2-ez2(i1))/dsqrt(z2var)))
      if(f2.le.fstar) then
         b2=-zrg
         f2=0.0
         do 220 i1=1,m1
220         f2=f2+h1(i1)*(1-norcdf((b2-ez2(i1))/dsqrt(z2var)))
         if(f2.lt.fstar) then
            zbdy(1,k2)=-zrg
            goto 280
         endif
      endif

c  Iteration to narrow the range (b1,b2).  Initially, we have
c
c     f1 < fstar < f2
c
c  and each iteration preserves this condition.

250   continue

c     print *,'250:b1,b2,f1,f2 ',b1,b2,f1,f2

      if(f2-f1 .le. epslon) goto 270

      b3=0.5*(b1+b2)
      f3=0.0
      do 260 i1=1,m1
260      f3=f3+h1(i1)*(1-norcdf((b3-ez2(i1))/dsqrt(z2var)))

      if(f3.lt.fstar) then
         b1=b3
         f1=f3
      else
         b2=b3
         f2=f3
      endif

      goto 250

270   zbdy(2,k2)=b1+(b2-b1)*(fstar-f1)/(f2-f1)

280   continue

      if(zbdy(1,k2).gt.zbdy(2,k2)) then
         write(6,290)
290      format(' Warning in subroutine ersp2:  Value calculated ',
     +          'for lower boundary',
     +         /'point is higher than value for upper boundary point.'
     +         /'The average of the two will be used for both.')
         write(6,291) zbdy(1,k2),zbdy(2,k2)
291      format(' zbdy(1,2;k2) =',2(2x,f10.6))
         zbdy(1,k2)=0.5*(zbdy(1,k2)+zbdy(2,k2))
         zbdy(2,k2)=zbdy(1,k2)
      endif

c  Update pu and pl

      do 400 i1=1,m1
         pu=pu+h1(i1)*(1-norcdf((zbdy(2,k2)-ez2(i1))/dsqrt(z2var)))
         pl=pl+h1(i1)*norcdf((zbdy(1,k2)-ez2(i1))/dsqrt(z2var))
400   continue
 
      if(k2.eq.na) goto 600

c   Calculate h2(i2) (=h_k2)

c   Set values z2(i2), i2=1,...,m2

c   Place an initial mesh of points suited to the N(theta*sq(inf2),1)
c   marginal distribution of Z2.

      zmean=0.0
      lim2=r-1
      do 410 i=1,lim2
410      zmesh(i)=zmean-(3.0+4.0*dlog(r/(1.0d0*i)))
      lim2=5*r
      do 420 i=r,lim2
420      zmesh(i)=zmean-(3.0-3.0d0*(i-r)/(2.0d0*r))
      lim1=5*r+1
      lim2=6*r-1
      do 430 i=lim1,lim2
430      zmesh(i)=zmean+(3.0+4.0*dlog(r/(6.0d0*r-i)))
      nmesh=6*r-1

c   Remove mesh points outside the range zbdy(1,k2) to zbdy(2,k2) and
c   add these points if appropriate.

      j1=1
      do 440 i=1,nmesh
         if(zmesh(i).lt.zbdy(1,k2)) j1=i
440      if(zmesh(i).lt.zbdy(1,k2) .and. i.lt.nmesh) zmesh(i)=zbdy(1,k2)
      j2=nmesh
      do 450 ii=1,nmesh
         i=nmesh+1-ii
         if(zmesh(i).gt.zbdy(2,k2)) j2=i
450      if(zmesh(i).gt.zbdy(2,k2) .and. i.gt.1) zmesh(i)=zbdy(2,k2)
      nmesh2=1+j2-j1

c   Create Simpson's rule grid points and associated weights.

      m2=2*nmesh2-1
      do 460 i=1,nmesh2
460      z2(2*i-1)=zmesh(j1+i-1)
      if(nmesh2.eq.1) goto 480
      lim2=nmesh2-1
      do 470 i=1,lim2
470      z2(2*i)=0.5*(z2(2*i-1)+z2(2*i+1))
480   continue
      w2(1)=0.0
      if(m2.eq.1) goto 520
      w2(1)=(z2(3)-z2(1))/6.0d0
      w2(m2)=(z2(m2)-z2(m2-2))/6.0d0
      if(m2.eq.3) goto 500
      lim2=m2-2
      do 490 i=3,lim2,2
490      w2(i)=(z2(i+2)-z2(i-2))/6.0d0
500   continue
      lim2=m2-1
      do 510 i=2,lim2,2
510      w2(i)=(z2(i+1)-z2(i-1))*4.0d0/6.0d0
520   continue

c   Calculate h2

      do 530 i2=1,m2
         h2(i2)=0.0
         do 540 i1=1,m1
540         h2(i2)=h2(i2)+f(z2(i2)-ez2(i1),z2var)*w2(i2)*h1(i1)
530      continue

c   Over-write k1,m1,z1 and h1 with new values

      k1=k2
      inf1=inf2
      m1=m2
      do 550 i=1,m1
         z1(i)=z2(i)
550      h1(i)=h2(i)

      goto 50 

600   continue

      goto 999

c  Error messages

980   write(6,981)
981   format(' Value of na out of range 1 to 100')
      goto 999

982   write(6,983)
983   format(' Negative value of inf(k)')
      goto 999

984   write(6,985)
985   format(' Decreasing values of inf(k)')
      goto 999

986   write(6,987) 
987   format(' delta =< 0')
      goto 999
986   write(6,987) 

988   write(6,989) 
989   format(' cume out of range or non-increasing')
      goto 999

990   write(6,991)
991   format(' alpha or beta =< 0 or >= 1')
      goto 999
 
992   write(6,993)
993   format(' ind not equal to 1 or 2')
      goto 999
 
994   write(6,995)
995   format(' Value of r out of range 1 to 16')
      goto 999
 
999   return
      end
