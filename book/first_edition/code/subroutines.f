c
c  Subroutines for normal responses
c

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

c=======================================================================
c
      subroutine gst2(na,inf,theta,zbdy,r,pu,pe,pl,pbya,einf,ierr)
c
c-----------------------------------------------------------------------
c
c   GST2 calculates properties of a group sequential two-sided test.
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
c   with zbdy(k) < zbdy(k) for k=1,...,na.
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
c      zbdy(1,k) < zbdy(2,k) for each k
c      r in the range 3 to 16.
c
c-----------------------------------------------------------------------
c
c   OUTPUT
c
c      pu         Pr{exit through upper boundary}               REAL*8
c      pe         Pr{remain within the boundary throughout}     REAL*8
c      pl         Pr{exit through lower boundary}               REAL*8
c
c      pbya(3,k)  Pr{Exit by upper boundary at analysis k}      REAL*8
c      pbya(2,k)  Pr{Terminate inside boundary at analysis k}   REAL*8
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
      real*8 inf(100),theta,zbdy(2,100),pu,pe,pl,pbya(3,100),einf,
     +       z1(200),z2(200),h1(200),h2(200),inf1,inf2,z2var,
     +       pupper,plower,ez2(200),zmean,zmesh(100),w2(200),
     +       norcdf,f,x,v

      f(x,v)=dexp(-x*x/(2.0d0*v))/dsqrt(6.2831853d0*v)

c   Check input variables

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

      do 4 k=1,na
4        if(zbdy(1,k).ge.zbdy(2,k)) ierr=1
      if(ierr.eq.1) goto 986

      if(r.lt.1 .or. r.gt.16) ierr=1
      if(ierr.eq.1) goto 988

c   Initialise probabilities and expected information

      pu=0.0
      pe=0.0
      pl=0.0
      do 10 k=1,na
      do 10 j=1,3
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
         pbya(3,k2)=pbya(3,k2)+pupper*h1(i1)
         if(k2.eq.na) then
            pe=pe+(1-pupper-plower)*h1(i1)
            pbya(2,k2)=pbya(2,k2)+(1-pupper-plower)*h1(i1)
         endif
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
989   format(' Value of r out of range 1 to 16')
      goto 999
 
999   return
      end
c=======================================================================
c
      subroutine ersp2(na,inf,cume,r,zbdy,ierr)
c
c-----------------------------------------------------------------------
c
c   ERSP2 computes a boundary for a group sequential two-sided test
c   according to the error spending approach. 
c
c                                            C. Jennison, 26 August 1999
c
c-----------------------------------------------------------------------
c
c   INPUT
c
c      Number of analyses = na                     INTEGER
c      Information levels = inf(1),...,inf(na)     REAL*8(100)
c
c   Cumulative error probabilities under theta=0 are specified as
c
c      exiting upper boundary = cume(2;1,...,na)
c      exiting lower boundary = cume(1;1,...,na)   REAL*8(2,100)
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
c   PERMISSIBLE  VALUES
c
c      0 < na =< 100
c      inf(k) > 0 for each k and inf(k) strictly increasing with k
c
c      cume(1,k) > 0 for each k and cume(1,k) strictly increasing with k
c      cume(2,k) > 0 for each k and cume(2,k) strictly increasing with k
c      cume(1,na) + cume(2,na) =< 1
c
c      r in the range 3 to 16.
c
c-----------------------------------------------------------------------
c
c   OUTPUT
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
      real*8 inf(100),cume(2,100),zbdy(2,100),
     +       z1(200),z2(200),h1(200),h2(200),inf1,inf2,ez2(200),z2var,
     +       pu,pl,b1,b2,b3,f1,f2,f3,fstar,zmean,zmesh(100),w2(200),
     +       zrg,epslon,
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

      do 3 k=1,na
         if(cume(2,k).le.0) ierr=1
         if(cume(1,k).le.0) ierr=1
3     continue
      if(ierr.eq.1) goto 986

      if(na.ge.2) then
         lim2=na-1
         do 4 k=1,lim2
            if(cume(2,k+1).le.cume(2,k)) ierr=1
            if(cume(1,k+1).le.cume(1,k)) ierr=1
4        continue
      endif
      if(ierr.eq.1) goto 986

      if((cume(1,na)+cume(2,na)).gt.1.0) ierr=1
      if(ierr.eq.1) goto 988

      if(r.lt.1 .or. r.gt.16) ierr=1
      if(ierr.eq.1) goto 990

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
987   format(' cume negative or non-increasing')
      goto 999

988   write(6,989)
989   format(' cume(1,na) + cume(2,na) > 1')
      goto 999

990   write(6,991)
991   format(' Value of r out of range 1 to 16')
      goto 999
 
999   return
      end

c=======================================================================
c
      real*8 function norcdf(x)
c
c-----------------------------------------------------------------------
c
c   NORCDF returns the cumulative distribution of a standard normal
c   distribution.
c
c   This function calls on the NAG library routine  S15ACF.
c   An equivalent routine must be substituted if you do not have access
c   to the NAG library.
c
c                                            C. Jennison, 25 August 1999
c
c-----------------------------------------------------------------------
c
c   INPUT     x                                     REAL*8
c
c   OUTPUT    norcdf    Pr{Z < x} when Z ~ N(0,1)   REAL*8
c
c-----------------------------------------------------------------------
c

      integer ifail
      real*8 x,s15acf

      ifail=0
      norcdf=1-s15acf(x,ifail)

      return
      end

c=======================================================================
c
      real*8 function invnor(p)
c
c-----------------------------------------------------------------------
c
c   INVNOR returns a quantile of the standard normal distribution.
c
c   This function calls on the NAG library routine  G01FAF.
c   An equivalent routine must be substituted if you do not have access
c   to the NAG library.
c
c                                            C. Jennison, 25 August 1999
c
c-----------------------------------------------------------------------
c
c   INPUT     p                                        REAL*8
c
c   OUTPUT    invnor    The value of x such that       REAL*8
c                       Pr{Z<x}=p when Z ~ N(0,1)
c
c-----------------------------------------------------------------------
c

      integer ifail
      real*8 p,g01faf

      ifail=0
      invnor=g01faf('L',p,ifail)

      return
      end

c=======================================================================
