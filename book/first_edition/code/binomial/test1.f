      integer n,nstar,ierr
      real*8 p,q,
     +       invnor,norcdf

      n=4
      p=0.5
c     q=1-norcdf(2.067d0)
      q=0.4
      print*,'q=',q

      call binv(n,p,q,nstar,ierr)

      print*,'n=',n
      print*,'p=',p
      print*,'q=',q
      print*,'nstar=',nstar

      stop
      end
