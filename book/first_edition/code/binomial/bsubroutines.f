c
c    Please notify Chris Jennison (cj@maths.bath.ac.uk) if you have any
c    problems in using these programs.
c
c                                           C. Jennison, 5 December 2001
c

c=======================================================================
c
      subroutine bgst(nobs,n,a,b,p,pu,pl,en,ierr)
c
c-----------------------------------------------------------------------
c
c   BGST calculates exit probabilities through upper and lower
c   boundaries for a group sequential test with Bernoulli responses,
c   each with success probability p.
c
c-----------------------------------------------------------------------
c
c   INPUT
c
c      Number of stages observed       = nobs    INTEGER
c      Group sizes           = n(1,...,nobs)     INTEGER(100)
c      Lower boundary        = a(1,...,nobs)     INTEGER(100)
c      Upper boundary        = b(1,...,nobs)     INTEGER(100)
c      Success prob. for a single obs. = p       REAL*8
c
c   At stage k the upper boundary is crossed if the number of
c   successes >= b(k), and the lower boudary is crossed if the number
c   of successes =< a(k).
c   Values a(k)=-1 and b(k)=n(1)+...+n(k)+1 denote that there is no
c   possibility of crossing the lower or upper boundary, respectively,
c   at stage k.
c   It is required that a(nobs)=b(nobs)-1 so that stopping by stage nobs
c   is certain.
c
c-----------------------------------------------------------------------
c
c   PERMISSIBLE  VALUES
c
c      0 < nobs =< 100
c      0 =< sobs =< 3000
c      0 =< p =< 1
c      n(1) > 0, ... , n(obs) > 0 and n(1)+...+n(nobs) =< 3000
c      -1 =< a(k) =< n(1)+...+n(k)
c       0 =< b(k) =< n(1)+...+n(k)+1
c      a(k) < b(k) for k < nobs,  a(nobs)=b(nobs)-1
c
c-----------------------------------------------------------------------
c
c   OUTPUT
c
c      pu      P(cross upper boundary)                          REAL*8
c      pl      P(cross lower boundary)                          REAL*8
c      en      E(number of observations at termination)         REAL*8
c      ierr    error indicator                                  INTEGER
c                    =1  if an error is detected
c                    =0  if no error occurs
c
c-----------------------------------------------------------------------
c
c   SUBROUTINE and FUNCTION calls
c
c      BPRCAL
c
c-----------------------------------------------------------------------
c

      real*8 p,pu,pl,en,c(3001),cc(3001),bp(3001),bple(3001),bpge(3001)
      integer nobs,n(100),a(100),b(100),ierr,m(100),k,n1,n2,ierr2,i

c   Check input

      ierr=0
      if(nobs.lt.1 .or. nobs.gt.100) ierr=1
      if(n(1).lt.1) ierr=1
      m(1)=n(1)
      if(nobs.eq.1) goto 20
      do 10 k=2,nobs
         if(n(k).lt.1) ierr=1
         m(k)=m(k-1)+n(k)
10       continue
20    if(m(nobs).gt.3000) ierr=1
      do 30 k=1,nobs
30       if(a(k).lt.-1 .or. b(k).gt.m(k)+1 .or. a(k).ge.b(k)) ierr=1
      if(a(nobs).ne.b(nobs)-1) ierr=1
      if(p.lt.0 .or. p.gt.1) ierr=1
      if(ierr.eq.1) goto 900

c   Initialise variables

      pu=0.0
      pl=0.0
      en=0.0

      k=0
      n1=1
      n2=1
      c(1)=1.0d0

100   continue

      call bprcal(n(k+1),p,bp,bple,bpge,ierr2)
      if(ierr2.eq.1) goto 910

      if(b(k+1).eq.m(k+1)+1) goto 115
      do 110 i=n1,n2
         if(i-1.lt.b(k+1)-n(k+1)) goto 110
         pu=pu+c(i)*bpge(1+max0(0,b(k+1)-i+1))
110      continue
115   continue

      if(a(k+1).eq.-1) goto 125
      do 120 i=n1,n2
         if(i-1.gt.a(k+1)) goto 120
         pl=pl+c(i)*bple(1+min0(a(k+1)-i+1,n(k+1)))
120      continue
125   continue

      do 130 i=n1,n2
130      en=en+c(i)*n(k+1)

      if(k+1.eq.nobs) goto 200

      nn1=a(k+1)+2
      nn2=b(k+1)

      do 150 j=nn1,nn2
         cc(j)=0
         do 140 i=n1,n2
            if(j-i.lt.0 .or. j-i.gt.n(k+1)) goto 140
            cc(j)=cc(j)+c(i)*bp(1+j-i)
140         continue
150      continue

      n1=nn1
      n2=nn2
      do 160 i=n1,n2
160      c(i)=cc(i)
      k=k+1

      goto 100

200   goto 999

900   write(6,901)
901   format(' Error in values entered to subroutine bgst')
      goto 999

910   write(6,911)
911   format(' Error in call to subroutine bprcal')
      ierr=2
      goto 999

999   return
      end

c=======================================================================
c
      subroutine bpci(nobs,sobs,n,a,b,p,pgte,plte,ierr)
c
c-----------------------------------------------------------------------
c
c   BPCI calculates the probabilities that the outcome of a group
c   sequential test with Bernoulli responses, each with success
c   probability p, should be less than or equal to or greater than or
c   equal to a particular observed outcome.
c
c                                            C. Jennison, 3 October 2000
c
c-----------------------------------------------------------------------
c
c   INPUT
c
c      Number of stages observed       = nobs    INTEGER
c      Number successes on termination = sobs    INTEGER
c      Success prob. for a single obs. = p       REAL*8
c      Group sizes           = n(1,...,nobs)     INTEGER(100)
c      Lower boundary        = a(1,...,nobs-1)   INTEGER(100)
c      Upper boundary        = b(1,...,nobs-1)   INTEGER(100)
c
c   At stage k the upper boundary is crossed if the number of
c   successes >= b(k), and the lower boudary is crossed if the number
c   of successes =< a(k).
c   Values a(k)=-1 and b(k)=n(1)+...+n(k)+1 denote that there is no
c   possibility of crossing the lower or upper boundary, respectively,
c   at stage k.
c
c-----------------------------------------------------------------------
c
c   PERMISSIBLE  VALUES
c
c      0 < nobs =< 100
c      0 =< sobs =< 3000
c      0 =< p =< 1
c      n(1) > 0, ... , n(obs) > 0 and n(1)+...+n(nobs) =< 3000
c      -1 =< a(k) =< n(1)+...+n(k)
c       0 =< b(k) =< n(1)+...+n(k)+1
c      a(k) < b(k) for k =< nobs
c
c-----------------------------------------------------------------------
c
c   OUTPUT
c
c      pgte    Pr{cross upper boundary or end with              REAL*8
c                  >= sobs successes at stage nobs}
c      plte    Pr{cross lower boundary or end with              REAL*8
c                  =< sobs successes at stage nobs}
c      ierr    error indicator                                  INTEGER
c                    =1  if an error is detected
c                    =0  if no error occurs
c
c-----------------------------------------------------------------------
c
c   SUBROUTINE and FUNCTION calls
c
c      BPRCAL
c
c-----------------------------------------------------------------------
c 

      real*8 p,pgte,plte,c(3001),cc(3001),bp(3001),bple(3001),
     +       bpge(3001),pu,pl
      integer nobs,sobs,n(100),a(100),b(100),ierr,m(100),i,j,k,n1,n2,
     +        ierr2,nn1,nn2,anext,bnext

c   Check input

      ierr=0
      if(nobs.lt.1 .or. nobs.gt.100) ierr=1
      if(n(1).lt.1) ierr=1
      m(1)=n(1)
      if(nobs.eq.1) goto 20
      do 10 k=2,nobs
         if(n(k).lt.1) ierr=1
         m(k)=m(k-1)+n(k)
10       continue
20    if(m(nobs).gt.3000) ierr=1
      do 30 k=1,nobs
30       if(a(k).lt.-1 .or. b(k).gt.m(k)+1 .or. a(k).ge.b(k)) ierr=1
      if(p.lt.0.0d0 .or. p.gt.1.0d0) ierr=1
      if(ierr.eq.1) goto 900

c   Initialise variables

      pu=0.0
      pl=0.0

      k=0
      n1=1
      n2=1
      c(1)=1.0d0

100   continue

      call bprcal(n(k+1),p,bp,bple,bpge,ierr2)

      if(ierr2.eq.1) goto 910

      if(k+1.lt.nobs) then
         anext=a(k+1)
         bnext=b(k+1)
      else
         anext=sobs
         bnext=sobs
      endif

      if(bnext.eq.m(k+1)+1) goto 115
      do 110 i=n1,n2
         if(i-1.lt.bnext-n(k+1)) goto 110
         pu=pu+c(i)*bpge(1+max0(0,bnext-i+1))
110      continue
115   continue

      if(k+1.lt.nobs .and. a(k+1).eq.-1) goto 125
      do 120 i=n1,n2
         if(i-1.gt.anext) goto 120
         pl=pl+c(i)*bple(1+min0(anext-i+1,n(k+1)))
120      continue
125   continue

      if(k+1.eq.nobs) goto 200

      nn1=a(k+1)+2
      nn2=b(k+1)

      do 150 j=nn1,nn2
         cc(j)=0
         do 140 i=n1,n2
            if(j-i.lt.0 .or. j-i.gt.n(k+1)) goto 140
            cc(j)=cc(j)+c(i)*bp(1+j-i)
140         continue
150      continue

      n1=nn1
      n2=nn2
      do 160 i=n1,n2
160      c(i)=cc(i)

      k=k+1

      goto 100

200   continue

      pgte=pu
      plte=pl

      goto 999

900   write(6,901)
901   format(' Error in values entered to subroutine bpci')
      goto 999

910   write(6,911)
911   format(' Error in call to subroutine bprcal')
      ierr=2
      goto 999

999   return
      end

c=======================================================================
c
      subroutine bprcal(n,p,bp,bple,bpge,ierr)
c
c-----------------------------------------------------------------------
c
c   BPRCAL calculates the probability that a Binomial(n,p) random
c   variable takes each individual value 0,...,n.  It also returns
c   cumulative probabilties of the variable being =< or >= each
c   possible value.
c
c                                            C. Jennison, 3 October 2000
c
c-----------------------------------------------------------------------
c
c   INPUT
c
c      Number of bernoulli trials      = n       INTEGER
c      Success prob. for a single obs. = p       REAL*8
c
c-----------------------------------------------------------------------
c
c   PERMISSIBLE  VALUES
c
c      0 < n =< 3000
c      0 =< p =< 1
c
c-----------------------------------------------------------------------
c
c   OUTPUT
c
c      bp(i=1,...,n+1)     = Pr{X = i-1}                    REAL*8(3001)
c      bple(i=1,...,n+1)   = Pr{X =< i-1}                   REAL*8(3001)
c      bpge(i=1,...,n+1)   = Pr{X >= i-1}                   REAL*8(3001)
c      ierr    error indicator                              INTEGER
c                    =1  if an error is detected
c                    =0  if no error occurs
c
c-----------------------------------------------------------------------
c

      real*8 p,bp(3001),bple(3001),bpge(3001),pmax,sum
      integer n,ierr,imax,i,j,lim

c  Check input

      ierr=0
      if(n.lt.0 .or. n.gt.3000) ierr=1
      if(p.lt.0.0d0 .or. p.gt.1.0d0) ierr=1
      if(ierr.eq.1) goto 900

      imax=min0(n,max0(0,ifix(sngl(n*p))))

      if(imax.eq.0) then
         pmax=(1.0d0-p)**n
         goto 11
      elseif(imax.eq.n) then
         pmax=p**n
         goto 11
      endif

      sum=0.0d0
      do 10 j=1,imax
10       sum=sum+dlog(n-j+1.0d0)-dlog(dfloat(j))
      sum=sum+imax*dlog(p)+(n-imax)*dlog(1.0d0-p)
      pmax=dexp(sum)

11    continue

      bp(imax+1)=pmax

      if(imax.eq.0) goto 21

c   find probabilities for 0, ... , imax-1

      do 20 j=1,imax
         bp(imax+1-j)=bp(imax+2-j)*(imax-j+1.0d0)*(1.0d0-p)
     +                /((n-imax+j)*p)
20       continue

21    continue

      if(imax.eq.n) goto 31

c   find probabilities of imax+1, ... , n

      lim=n-imax
      do 30 j=1,lim
         bp(imax+1+j)=bp(imax+j)*(n-imax-j+1.0d0)*p
     +                /((imax+j)*(1.0d0-p))
30       continue

31    continue

      bple(1)=bp(1)
      bpge(1)=1.0d0
      lim=n+1
      do 40 i=2,lim
         bple(i)=bple(i-1)+bp(i)
40       bpge(i)=bpge(i-1)-bp(i-1)

      goto 999

900   write(6,901) n,p
901   format(' Error in subroutine bprcal')
      goto 999

999   return
      end

c=======================================================================
c
      subroutine findci(ns,n,a,b,conf,nobs,sobs,ci,ierr)
c
c-----------------------------------------------------------------------
c
c   FINDCI calculates a confidence interval for a success probability,
c   p, on the termination of a group sequential test with Bernoulli
c   responses.
c
c                                            C. Jennison, 3 October 2000
c
c-----------------------------------------------------------------------
c
c   INPUT
c
c      Number of stages                    = ns      INTEGER
c      Stage at which termination occured  = nobs    INTEGER
c      Number successes on termination     = sobs    INTEGER
c      Group sizes           = n(1,...,ns)           INTEGER(100)
c      Lower boundary        = a(1,...,ns)           INTEGER(100)
c      Upper boundary        = b(1,...,ns)           INTEGER(100)
c      Confidence level of the interval    = conf    REAL*8
c
c   At stage k the upper boundary is crossed if the number of
c   successes >= b(k), and the lower boudary is crossed if the number
c   of successes =< a(k).
c   Values a(k)=-1 and b(k)=n(1)+...+n(k)+1 denote that there is no
c   possibility of crossing the lower or upper boundary, respectively,
c   at stage k.
c
c   The end points of the desired confidence interval are:
c
c      lower endpoint ci(1) = p such that
c
c        Pr{Obtain observed outcome or higher} = (1-conf)/2
c
c      upper endpoint ci(2) = p such that
c
c        Pr{Obtain observed outcome or lower} = (1-conf)/2
c
c   Here, the ordering is defined by, in order of priority,
c
c      (1) the boundary crossed,
c      (2) the stage reached,
c      (3) the number of successes observed.
c
c-----------------------------------------------------------------------
c
c   PERMISSIBLE  VALUES
c
c      0 < ns =< 100
c      0 < nobs =< 100
c      0 =< sobs =< 3000
c      n(1) > 0, ... , n(obs) > 0 and n(1)+...+n(ns) =< 3000
c      -1 =< a(k) =< n(1)+...+n(k)
c       0 =< b(k) =< n(1)+...+n(k)+1
c      a(k) < b(k) for k < ns and a(ns)=b(ns)-1
c      0 < conf < 1
c
c-----------------------------------------------------------------------
c
c   OUTPUT
c
c     ci(1)   Lower end point of confidence interval            REAL*8
c     ci(2)   Upper end point of confidence interval            REAL*8
c     ierr    error indicator                                   INTEGER
c                   =0 if all is OK,
c                   =1 if input values out of range,
c                   =2 if error in call to subroutine bpci
c
c-----------------------------------------------------------------------
c
c   SUBROUTINE and FUNCTION calls
c
c      BPCI
c
c-----------------------------------------------------------------------
c 

      integer ns,n(100),a(100),b(100),nobs,sobs,ierr,m(100),k,ierr2,
     +        ns1low,ns1upp,ilow,iupp
      real*8 conf,ci(2),p1,p2,pnew,f1,f2,fnew,pgte,plte

c   Check input values

      ierr=0
      if(ns.lt.1 .or. ns.gt.100) ierr=1
      if(n(1).lt.1) ierr=1
      m(1)=n(1)
      if(ns.eq.1) goto 20
      do 10 k=2,ns
         if(n(k).lt.1) ierr=1
         m(k)=m(k-1)+n(k)
10       continue
20    if(m(ns).gt.3000) ierr=1
      do 30 k=1,ns
30       if(a(k).lt.-1 .or. b(k).gt.m(k)+1 .or. a(k).ge.b(k)) ierr=1
      if(a(ns).ne.b(ns)-1) ierr=1
      if(sobs.lt.b(nobs) .and. sobs.gt.a(nobs)) ierr=1
      if(ierr.eq.1) goto 900

      ilow=0
      iupp=0
      do 40 k=1,ns
         if(ilow.eq.0 .and. a(k).gt.-1) then
            ns1low=k
            ilow=1
         endif
         if(iupp.eq.0 .and. b(k).lt.m(k)+1) then
            ns1upp=k
            iupp=1
         endif
40    continue

c  Finding ci(1)

      if(nobs.eq.ns1low .and. sobs.eq.0) then
         ci(1)=0.0
         goto 110
      endif

c  Initialise bisection search for root of pu = (1-conf)/2

      p1=0.0
      call bpci(nobs,sobs,n,a,b,p1,pgte,plte,ierr2)
      if(ierr.eq.2) goto 910
      f1=pgte-(1-conf)/2

      p2=1.0
      call bpci(nobs,sobs,n,a,b,p2,pgte,plte,ierr2)
      if(ierr.eq.2) goto 910
      f2=pgte-(1-conf)/2

c  bisection search

50    continue

      if(dabs(f1-f2).le.0.00001) goto 100
      pnew=(p1+p2)/2
      call bpci(nobs,sobs,n,a,b,pnew,pgte,plte,ierr2)
      if(ierr.eq.2) goto 910
      fnew=pgte-(1-conf)/2

      if(fnew.ge.0) then
         p2=pnew
         f2=fnew
      else
         p1=pnew
         f1=fnew
      endif

      goto 50

100   continue

      ci(1)=(p1*f2-p2*f1)/(f2-f1)

110   continue

c  Finding ci(2)

c  Initialise binary search for root of pl = (1-conf)/2

      if(nobs.eq.ns1upp .and. sobs.eq.m(nobs)) then
         ci(2)=1.0
         goto 210
      endif

      p1=0.0
      call bpci(nobs,sobs,n,a,b,p1,pgte,plte,ierr2)
      if(ierr.eq.2) goto 910
      f1=plte-(1-conf)/2

      p2=1.0
      call bpci(nobs,sobs,n,a,b,p2,pgte,plte,ierr2)
      if(ierr.eq.2) goto 910
      f2=plte-(1-conf)/2

150   continue

      if(dabs(f1-f2).le.0.00001) goto 200
      pnew=(p1+p2)/2
      call bpci(nobs,sobs,n,a,b,pnew,pgte,plte,ierr2)
      if(ierr.eq.2) goto 910
      fnew=plte-(1-conf)/2

      if(fnew.le.0) then
         p2=pnew
         f2=fnew
      else
         p1=pnew
         f1=fnew
      endif

      goto 150

200   continue

      ci(2)=(p1*f2-p2*f1)/(f2-f1)

210   continue

      goto 999

900   write(6,901)
901   format(' Error in values entered to subroutine findci')
      goto 999

910   write(6,911)
911   format(' Error in call to subroutine bpci')
      ierr=2
      goto 999

999   return
      end
