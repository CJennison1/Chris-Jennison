      integer na,r,ierr,k
      real*8 inf(100),xvar,cume(2,100),sbdy(2,100),zbdy(2,100),
     +       alf(2,100)

      na=5
      r=16

      xvar=1.0

      do 5 k=1,na
         inf(k)=k
         cume(2,k)=(k*0.1)/na
         cume(1,k)=(k*0.05)/na
5     continue

      write(6,6) (cume(2,k),k=1,na)
6     format('Upper error to spend ',10f8.4)
      write(6,7) (cume(1,k),k=1,na)
7     format('Lower error to spend ',10f8.4)

      call ldtest(na,inf,xvar,cume,r,sbdy,alf,ierr)

      write(6,11)ierr
11    format(/'  ierr = ',i2)
      write(6,12) (sbdy(2,k),k=1,na)
12    format(/'  Upper s boundary ',10f8.4)
      write(6,13) (sbdy(1,k),k=1,na)
13    format('  Lower s boundary ',10f8.4)
      write(6,14) (alf(2,k),k=1,na)
14    format(/'  Upper error spent ',10f10.6)
      write(6,15) (alf(1,k),k=1,na)
15    format('  Lower error spent ',10f10.6)

      print *,' '

      write(6,6) (cume(2,k),k=1,na)
      write(6,7) (cume(1,k),k=1,na)

      call ersp2(na,inf,cume,r,zbdy,ierr)

      do 20 k=1,na
         sbdy(2,k)=dsqrt(inf(k))*zbdy(2,k)
         sbdy(1,k)=dsqrt(inf(k))*zbdy(1,k)
20    continue

      write(6,11)ierr
      write(6,12) (sbdy(2,k),k=1,na)
      write(6,13) (sbdy(1,k),k=1,na)
      write(6,16) (zbdy(2,k),k=1,na)
16    format(/'  Upper z bdy ',10f8.4)
      write(6,17) (zbdy(1,k),k=1,na)
17    format('  Lower z bdy ',10f8.4)

      stop
      end
