      integer ns,n(100),a(100),b(100),sobs,ierr,k,m(100),i,
     +        lim1,lim2,table1(100,3001)
      real*8 conf,ci(2),table2(100,3001,2)

c     Example: MIL STD 105D, plan Q 1.0

      data ns / 7 /
      data n  / 7*315, 93*0 /
      data b  / 9, 14, 19, 25, 29, 33, 38, 93*0 /
      data a  / 2,  7, 13, 19, 25, 31, 37, 93*0 /

      data conf / 0.90 /

      open(6,file='cis.out_1',fileopt='eof')

      do 1 k=1,ns
         if(k.eq.1) m(1)=n(1)
         if(k.gt.1) m(k)=m(k-1)+n(k)
1     continue

      do 2 k=1,ns
      do 2 i=0,m(k)
2        table1(k,i+1)=0

      do 100 k=1,ns

         if(a(k).gt.-1) then
            if(k.eq.1) then
               lim1=0
            else
               lim1=max0(a(k-1)+1,0)
            endif
            lim2=a(k)
            lim1=max0(lim1,lim2-10)
            if(lim1.gt.lim2) goto 20
            do 10 sobs=lim1,lim2
               call findci(ns,n,a,b,conf,k,sobs,ci,ierr)
               table1(k,sobs+1)=1
               table2(k,sobs+1,1)=ci(1)*100
               table2(k,sobs+1,2)=ci(2)*100
10          continue
20          continue
         endif

         if(b(k).lt.m(k)+1) then
            lim1=b(k)
            lim2=m(k)
            lim2=min0(lim2,lim1+10)
            do 50 sobs=lim1,lim2
               call findci(ns,n,a,b,conf,k,sobs,ci,ierr)
               table1(k,sobs+1)=1
               table2(k,sobs+1,1)=ci(1)*100
               table2(k,sobs+1,2)=ci(2)*100
50          continue
         endif

100   continue

      do 200 k=1,ns
      do 200 sobs=0,m(k)
         if(table1(k,sobs+1).eq.1) then
            write(6,210) k,sobs,table2(k,sobs+1,1),table2(k,sobs+1,2)
210         format('stage ',i2,' no. successes ',i4,' ci ',f5.1,1x,f5.1)
         endif
200   continue

      stop
      end
