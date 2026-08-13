      integer ns,n(100),m(100),a(100),b(100),i,k,ierr,nstar
      real*8 qq,pval(2),p,pu,pl,en

c     Example: 4 stage test of p=0.05 vs p > 0.05 with significnace
c              levels from a Pocock test for one-sided Type I error
c              rate 0.05.

      data ns / 4 /
      data n  / 25, 25, 25, 25, 96*0 /

      data qq / 0.0192d0 /

c     The one-sided significance level for a 4-group Pocock test
c     with alpha=0.05, C_p=2.07.

      data pval / 0.05d0, 0.1d0 /

c     open(6,file='example_4.out',fileopt='eof')

      print *,'qq=',qq

      p=pval(1)

      do 10 k=1,ns

         if(k.eq.1) then
            m(k)=n(1)
         else
            m(k)=m(k-1)+n(k)
         endif

         call binv(m(k),p,qq,nstar,ierr)

         b(k)=nstar
         if(k.ne.ns) then
            a(k)=-1
         else
            a(k)=b(k)-1
         endif

10    continue

         write(6,14)
14       format('Stage   Boundary points'/)
      do 20 k=1,ns
         write(6,15) k,a(k),b(k)
15       format(i3,2x, 3x, 2x,i3,',',2x,i3)
20    continue

      do 50 i=1,2

      p=pval(i)
      call bgst(ns,n,a,b,p,pu,pl,en,ierr)

      if(ierr.eq.1) write(6,30)
30    format('Error in bgst')
      if(ierr.eq.1) goto 40
      write(6,35) p,pu,pl,en
35    format(/' p = ',f6.4,5x,' p(exit upper boundary) = ',f8.5
     +       /16x,' p(exit lower boundary) = ',f8.5
     +       /16x,'         E(sample size) = ',f8.4)
40    continue
 
50    continue

      stop
      end
