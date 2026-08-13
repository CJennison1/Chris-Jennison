      integer ns,n(100),a(100),b(100),i,ierr
      real*8 pval(2),p,pu,pl,en

c     Example: Fleming (1982) 3 stage test of p=0.1 vs p=0.3

      data ns / 3 /
      data n  / 15, 10, 10, 97*0 /
      data b  /  5, 6, 7, 97*0 /
      data a  /  0, 3, 6, 97*0 /

      data pval / 0.1d0, 0.3d0 /

      open(6,file='props.out_1',fileopt='eof')

      do 50 i=1,2

      p=pval(i)
      call bgst(ns,n,a,b,p,pu,pl,en,ierr)

      if(ierr.eq.1) write(6,10)
10    format('Error in bgst')
      if(ierr.eq.1) goto 30
      write(6,20) p,pu,pl,en
20    format(/' p = ',f6.4,5x,' p(exit upper boundary) = ',f8.5
     +       /16x,' p(exit lower boundary) = ',f8.5
     +       /16x,'         E(sample size) = ',f8.4)
30    continue
 
50    continue

      stop
      end
