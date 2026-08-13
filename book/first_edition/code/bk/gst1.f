
c=======================================================================
c
      subroutine gst1(na,inf,theta,zbdy,r,pu,pl,pbya,einf,ierr)
c
c-----------------------------------------------------------------------
c
c   GST1 calculates properties of a group sequential one-sided test.
c
c                                            C. Jennison, 25 August 1999
c
c-----------------------------------------------------------------------
c
c   INPUT
c
c      Number of analyses = na                   INTEGER
c      Information levels = inf(1),...,inf(na)   REAL*8(100)
c      Parameter value    = theta                REAL*8
c
c   Standardized statistics Z(1),...,Z(na) are assumed to follow the
c   canonical joint distribution
c
c      (Z(1),...,Z(na)) ~ multivariate normal,
c      Z(k) ~ N(theta sq(inf(k)),1),   k=1,...,na
c      cov(Z(k1),Z(k2)) = sq(inf(k1)/inf(k2)),   for k1 =< k2.
c
c   The stopping boundary is specified as
c
c      upper:  zbdy(2,1),...,zbdy(2,na)          REAL*8(100)
c      lower:  zbdy(1,1),...,zbdy(1,na)          REAL*8(100)
c
c   with zbdy(1,k) < zbdy(2,k) for k=1,...,na-1 and
c   zbdy(1,na) = zbdy(2,na).
c
c   The parameter 
c
c      r                                         INTEGER
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
c   PERMISSIBLE  VALUES
c
c      0 < na =< 100
c      inf(k) > 0 for each k and inf(k) strictly increasing with k
c      zbdy(1,k) < zbdy(2,k) for each k=1,...,na-1 and
c      zbdy(1,na) = zbdy(2,na), at least within 10**-5.
c      r in the range 3 to 16.
c
c-----------------------------------------------------------------------
c
c   OUTPUT
c
c      pu         Pr{exit through upper boundary}               REAL*8
c      pl         Pr{exit through lower boundary}               REAL*8
c
c      pbya(2,k)  Pr{Exit by upper boundary at analysis k}      REAL*8
c      pbya(1,k)  Pr{Exit by lower boundary at analysis k}      REAL*8
c
c      einf       E(inf) on termination                         REAL*8
c
c      ierr       error indicator                               INTEGER
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
      real*8 inf(100),theta,zbdy(2,100),pu,pl,pbya(2,100),einf,
     +       z1(200),z2(200),h1(200),h2(200),inf1,inf2,z2var,
     +       pupper,plower,ez2(200),zmean,zmesh(100),w2(200),
     +       norcdf,f,x,v

      f(x,v)=dexp(-x*x/(2.0d0*v))/dsqrt(6.2831853d0*v)

c   Check input variables.  Set final boundary points to be equal if
c   they are slightly different (e.g., beacause of some rounding error).

      ierr=0

      if(na.lt.1 .or. na.gt.100) ierr=1
      if(ierr.eq.1) goto 980

      do 1 k=1,na
1        if(inf(k).le.0) ierr=1
      if(ierr.eq.1) goto 982

      if(na.eq.1) goto 3
      do 2 k=2,na
2        if(inf(k).le.inf(k-1)) ierr=1
      if(ierr.eq.1) goto 984
3     continue

      if(na.ge.2) then
         lim2=na-1
         do 4 k=1,lim2
4           if(zbdy(1,k).ge.zbdy(2,k)) ierr=1
      endif
      if(ierr.eq.1) goto 986
      if(dabs(zbdy(1,na)-zbdy(2,na)).gt.1.0d-5) ierr=1
      if(ierr.eq.1) goto 988
      zbdy(1,na)=zbdy(2,na)

      if(r.lt.1 .or. r.gt.16) ierr=1
      if(ierr.eq.1) goto 990

c   Initialise probabilities and expected information

      pu=0.0
      pl=0.0
      do 10 k=1,na
      do 10 j=1,2
         pbya(j,k)=0.0
10    continue
      einf=0.0

c   Iterative calculations
 
c   An integral at analysis k of e(z_k) over the sub-distribution of z_k,
c   for paths continuing to analysis k, is approximated by the sum
c
c      sum(i=1,m_k)  h_k(i_k) e(z_k(i_k))
c
c   The general step is to compute the grid of z values and associated
c   weights
c
c      z_k(i) and h_k(i),   i=1,m_k
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

c   Loop step from k1 to k2=k1+1

50    continue

      k2=k1+1
      inf2=inf(k2)

      z2var=1.0-inf1/inf2

c   Calculate ez2(i1), the conditional mean of Z2 given Z1 = z1(i1),
c   and update exit probabilities.

      do 100 i1=1,m1
         ez2(i1)=z1(i1)*dsqrt(inf1/inf2)+
     +           theta*((inf2-inf1)/dsqrt(inf2))
         pupper=1-norcdf((zbdy(2,k2)-ez2(i1))/dsqrt(z2var))
         plower=  norcdf((zbdy(1,k2)-ez2(i1))/dsqrt(z2var))
         pu=pu+pupper*h1(i1)
         pbya(2,k2)=pbya(2,k2)+pupper*h1(i1)
         pl=pl+plower*h1(i1)
         pbya(1,k2)=pbya(1,k2)+plower*h1(i1)
         einf=einf+(inf2-inf1)*h1(i1)
100      continue
 
      if(k2.eq.na) goto 300

c   Calculate h2(i2) (=h_k2)

c   Set values z2(i2), i2=1,...,m2

c   Place an initial mesh of points suited to the N(theta*sq(inf2),1)
c   marginal distribution of Z2.

      zmean=theta*dsqrt(inf2)
      lim2=r-1
      do 110 i=1,lim2
110      zmesh(i)=zmean-(3.0+4.0*dlog(r/(1.0d0*i)))
      lim2=5*r
      do 120 i=r,lim2
120      zmesh(i)=zmean-(3.0-3.0d0*(i-r)/(2.0d0*r))
      lim1=5*r+1
      lim2=6*r-1
      do 130 i=lim1,lim2
130      zmesh(i)=zmean+(3.0+4.0*dlog(r/(6.0d0*r-i)))
      nmesh=6*r-1

c   Remove mesh points outside the range zbdy(1,k2) to zbdy(2,k2) and
c   add these points if appropriate.

      j1=1
      do 140 i=1,nmesh
         if(zmesh(i).lt.zbdy(1,k2)) j1=i
140      if(zmesh(i).lt.zbdy(1,k2) .and. i.lt.nmesh) zmesh(i)=zbdy(1,k2)
      j2=nmesh
      do 150 ii=1,nmesh
         i=nmesh+1-ii
         if(zmesh(i).gt.zbdy(2,k2)) j2=i
150      if(zmesh(i).gt.zbdy(2,k2) .and. i.gt.1) zmesh(i)=zbdy(2,k2)
      nmesh2=1+j2-j1

c   Create Simpson's rule grid points and associated weights.

      m2=2*nmesh2-1
      do 160 i=1,nmesh2
160      z2(2*i-1)=zmesh(j1+i-1)
      if(nmesh2.eq.1) goto 180
      lim2=nmesh2-1
      do 170 i=1,lim2
170      z2(2*i)=0.5*(z2(2*i-1)+z2(2*i+1))
180   continue
      w2(1)=0.0
      if(m2.eq.1) goto 220
      w2(1)=(z2(3)-z2(1))/6.0d0
      w2(m2)=(z2(m2)-z2(m2-2))/6.0d0
      if(m2.eq.3) goto 200
      lim2=m2-2
      do 190 i=3,lim2,2
190      w2(i)=(z2(i+2)-z2(i-2))/6.0d0
200   continue
      lim2=m2-1
      do 210 i=2,lim2,2
210      w2(i)=(z2(i+1)-z2(i-1))*4.0d0/6.0d0
220   continue

c   Calculate h2

      do 230 i2=1,m2
         h2(i2)=0.0
         do 240 i1=1,m1
240         h2(i2)=h2(i2)+f(z2(i2)-ez2(i1),z2var)*w2(i2)*h1(i1)
230      continue

c   Over-write k1,m1,z1 and h1 with new values

      k1=k2
      inf1=inf2
      m1=m2
      do 250 i=1,m1
         z1(i)=z2(i)
250      h1(i)=h2(i)

      goto 50 

300   continue
      goto 999

c   Error messages

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
987   format(' Lower boundary zbdy(1,k) > upper boundary zbdy(2,k)')
      goto 999

988   write(6,989) 
989   format(' Final boundary points zbdy(1,na) and zbdy(2,na)'
     +       ' are not equal',
     +      /' (in fact they are not within 10**-5 of each other.')
      goto 999

990   write(6,991)
991   format(' Value of r out of range 1 to 16')
      goto 999
 
999   return
      end
