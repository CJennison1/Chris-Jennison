c=======================================================================
c
      subroutine binv(n,p,q,nstar,ierr)
c
c-----------------------------------------------------------------------
c
c   BINV finds the least value nstar such that the probability that a
c   Binomial(n,p) random variable >= nstar is less than or equal to q.
c
c                                           C. Jennison, 5 December 2001
c
c-----------------------------------------------------------------------
c
c   INPUT
c
c      Number of bernoulli trials      = n       INTEGER
c      Success prob. for a single obs. = p       REAL*8
c      Upper tail probability          = q       REAL*8
c
c-----------------------------------------------------------------------
c
c   PERMISSIBLE  VALUES
c
c      0 < n =< 3000
c      0 =< p =< 1
c      0 =< q =< 1
c
c-----------------------------------------------------------------------
c
c   OUTPUT
c
c      nstar   = value such that                               INTEGER
c
c                         p(Bin(n,p) >= nstar-1) > q 
c                and
c                         p(Bin(n,p) >= nstar) =< q,
c
c                which means nstar is set as n+1 if 
c
c                         p(Bin(n,p) = n) > q 
c
c      ierr    error indicator                                 INTEGER
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

      real*8 p,q,bp(3001),bple(3001),bpge(3001)
      integer n,nstar,ierr,ierr2,i

c  Check input

      ierr=0
      if(n.lt.0 .or. n.gt.3000) ierr=1
      if(p.lt.0.0d0 .or. p.gt.1.0d0) ierr=1
      if(q.lt.0.0d0 .or. q.gt.1.0d0) ierr=1
      if(ierr.eq.1) goto 900

      call bprcal(n,p,bp,bple,bpge,ierr2)
      if(ierr2.eq.1) goto 910

c       do 5 i=0,n
c          print *,'i=',i,' p(bin(n,p) >= i)=',bpge(i+1)
c 5     continue

      nstar=0
      do 10 i=0,n
         if(bpge(i+1).gt.q) nstar=i+1
10    continue
      goto 999

900   write(6,901)
901   format(' Error in values entered to subroutine binv')
      goto 999

910   write(6,911)
911   format(' Error in call to subroutine bprcal')
      ierr=2
      goto 999

999   return
      end
