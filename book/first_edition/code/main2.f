c-----------------------------------------------------------------------
c
c   MAIN2  Main program for calculating properties of selected two-sided
c   group sequential tests.
c
c                                            C. Jennison, 25 August 1999
c
c-----------------------------------------------------------------------
c
c   Three examples are presented
c
c     p)  a Pocock test
c     b)  an O'Brien & Fleming test
c     w)  a Wang & Tsiatis test with parameter Del=0.25
c
c   Constants defining these tests are taken from Tables
c
c         2.1 and 2.2 for the Pocock test
c         2.3 and 2.4 for the O'Brien & Fleming test
c         2.9 and 2.10 for the Wang & Tsiatis test
c
c   in the book Group Sequential Methods with Applications to Clinical
c   Trials by Jennison & Turnbull.
c
c   The expected sample sizes of these tests agree with the values
c   displayed in Tables 2.7, 2.8 and 2.12.  (Error probabilities do
c   differ a little from specified values as constants c and r are
c   rounded to 3 decimal places; similarly, E{inf} may differ in the
c   last decimal place from the tabulated value, which is based on
c   slightly more accurate values of c and r.)
c
c   For the Pocock and O'Brien & Fleming tests, probabilities of
c   stopping at each stage are similar to the probabilities shown in
c   Table 2.6.  Small differences arise because, there, information
c   levels were increased in rounding up to integer group sizes.   
c
c-----------------------------------------------------------------------

      integer na,r,ierr,k
      real*8 inf(100),theta,zbdy(2,100),pu,pe,pl,pbya(3,100),einf,
     +       alpha,beta,delta,ifix,cp,rp,cb,rb,cw,rw,delw,eratio,t1err,
     +       invnor

c   Set number of analyses, Type I and II error probabilities, and value
c   of theta at which power 1-beta is required.

      na=5
      alpha=0.05
      beta=0.1
      delta=2.0

c   Set the parameter controlling density of grid points in numerical
c   integrations.

      r=16

c   Set constants for the tests.

c   cp, cb and cw define the boundaries for the threee tests.
c   rp, rb and rw are ratios of the maximum information level in each
c   test to that needed in a fixed sample test.

      cp=2.413
      rp=1.207

      cb=2.040
      rb=1.026

      cw=2.136
      rw=1.066
      delw=0.25

c   Find information level needed by a fixed sample test

      ifix=(invnor(1-alpha/2)+invnor(1-beta))**2/delta**2

      write(6,1) alpha,beta,delta,ifix
1     format(/'A fixed sample, one-sided test with Type I'
     +        ' error probability =',f6.3,
     +       /'and Type II error probability =',f6.3,' at theta =',f6.3,
     +        ' requires an'
     +       /'information level =',f6.3,'.')

c   p) Pocock test

      do 110 k=1,na
         inf(k)=(k*rp*ifix)/na
         zbdy(2,k)= cp
         zbdy(1,k)=-zbdy(2,k)
110    continue

c   Under the null hypothesis

      theta=0.0
      call gst2(na,inf,theta,zbdy,r,pu,pe,pl,pbya,einf,ierr)
      eratio=einf/ifix
      t1err=pu+pl
      write(6,121) ierr,pu,pe,pl,t1err,einf,eratio
121   format(//'Pocock test'
     +        /'-----------'
     +       //'theta=0:'
     +       //'ierr=',i1
     +       //'Pr{Exit upper boundary}       = ',f10.5
     +        /'Pr{Continue inside boundary}  = ',f10.5
     +        /'Pr{Exit lower boundary}       = ',f10.5
     +       //'Type I error probability      = ',f10.5
     +       //'E{Information on termination} = ',f10.5
     +        /'E{Inf}/Fixed sample inf level = ',f10.5)
      write(6,122)
122   format(/'Probability of stopping by analysis'/)
      write(6,123) (pbya(3,k),k=1,na)
123   format('  upper boundary     ',5f10.6)
      write(6,124) (pbya(2,k),k=1,na)
124   format('  between boundaries ',5f10.6)
      write(6,125) (pbya(1,k),k=1,na)
125   format('  lower boundary     ',5f10.6)

c   At the positive alternative where power 1-beta is required

      theta=delta
      call gst2(na,inf,theta,zbdy,r,pu,pe,pl,pbya,einf,ierr)
      eratio=einf/ifix
      write(6,131) ierr,pu,pe,pl,pu,einf,eratio
131   format( /'theta=delta:'
     +       //'ierr=',i1
     +       //'Pr{Exit upper boundary}       = ',f10.5
     +        /'Pr{Continue inside boundary}  = ',f10.5
     +        /'Pr{Exit lower boundary}       = ',f10.5
     +       //'Power at theta=delta          = ',f10.5
     +       //'E{Information on termination} = ',f10.5
     +        /'E{Inf}/Fixed sample inf level = ',f10.5)
      write(6,132)
132   format(/'Probability of stopping by analysis'/)
      write(6,133) (pbya(3,k),k=1,na)
133   format('  upper boundary     ',5f10.6)
      write(6,134) (pbya(2,k),k=1,na)
134   format('  between boundaries ',5f10.6)
      write(6,135) (pbya(1,k),k=1,na)
135   format('  lower boundary     ',5f10.6)

c   b) O'Brien & Fleming test

      do 210 k=1,na
         inf(k)=(k*rb*ifix)/na
         zbdy(2,k)= cb*dsqrt(na*1.0d0/k)
         zbdy(1,k)=-zbdy(2,k)
210    continue

c   Under the null hypothesis

      theta=0.0
      call gst2(na,inf,theta,zbdy,r,pu,pe,pl,pbya,einf,ierr)
      eratio=einf/ifix
      t1err=pu+pl
      write(6,221) ierr,pu,pe,pl,t1err,einf,eratio
221   format(//'O''Brien & Fleming test'
     +        /'----------------------'
     +       //'theta=0:'
     +       //'ierr=',i1
     +       //'Pr{Exit upper boundary}       = ',f10.5
     +        /'Pr{Continue inside boundary}  = ',f10.5
     +        /'Pr{Exit lower boundary}       = ',f10.5
     +       //'Type I error probability      = ',f10.5
     +       //'E{Information on termination} = ',f10.5
     +        /'E{Inf}/Fixed sample inf level = ',f10.5)
      write(6,222)
222   format(/'Probability of stopping by analysis'/)
      write(6,223) (pbya(3,k),k=1,na)
223   format('  upper boundary     ',5f10.6)
      write(6,224) (pbya(2,k),k=1,na)
224   format('  between boundaries ',5f10.6)
      write(6,225) (pbya(1,k),k=1,na)
225   format('  lower boundary     ',5f10.6)

c   At the positive alternative where power 1-beta is required

      theta=delta
      call gst2(na,inf,theta,zbdy,r,pu,pe,pl,pbya,einf,ierr)
      eratio=einf/ifix
      write(6,231) ierr,pu,pe,pl,pu,einf,eratio
231   format( /'theta=delta:'
     +       //'ierr=',i1
     +       //'Pr{Exit upper boundary}       = ',f10.5
     +        /'Pr{Continue inside boundary}  = ',f10.5
     +        /'Pr{Exit lower boundary}       = ',f10.5
     +       //'Power at theta=delta          = ',f10.5
     +       //'E{Information on termination} = ',f10.5
     +        /'E{Inf}/Fixed sample inf level = ',f10.5)
      write(6,232)
232   format(/'Probability of stopping by analysis'/)
      write(6,233) (pbya(3,k),k=1,na)
233   format('  upper boundary     ',5f10.6)
      write(6,234) (pbya(2,k),k=1,na)
234   format('  between boundaries ',5f10.6)
      write(6,235) (pbya(1,k),k=1,na)
235   format('  lower boundary     ',5f10.6)

c   w) Wang & Tsiatis (Del=0.25) test

      do 310 k=1,na
         inf(k)=(k*rw*ifix)/na
         zbdy(2,k)= cw*(k*1.0d0/na)**(delw-0.5)
         zbdy(1,k)=-zbdy(2,k)
310    continue

c   Under the null hypothesis

      theta=0.0
      call gst2(na,inf,theta,zbdy,r,pu,pe,pl,pbya,einf,ierr)
      eratio=einf/ifix
      t1err=pu+pl
      write(6,321) ierr,pu,pe,pl,t1err,einf,eratio
321   format(//'Wang & Tsiatis (Del=0.25) test'
     +        /'------------------------------'
     +       //'theta=0:'
     +       //'ierr=',i1
     +       //'Pr{Exit upper boundary}       = ',f10.5
     +        /'Pr{Continue inside boundary}  = ',f10.5
     +        /'Pr{Exit lower boundary}       = ',f10.5
     +       //'Type I error probability      = ',f10.5
     +       //'E{Information on termination} = ',f10.5
     +        /'E{Inf}/Fixed sample inf level = ',f10.5)
      write(6,322)
322   format(/'Probability of stopping by analysis'/)
      write(6,323) (pbya(3,k),k=1,na)
323   format('  upper boundary     ',5f10.6)
      write(6,324) (pbya(2,k),k=1,na)
324   format('  between boundaries ',5f10.6)
      write(6,325) (pbya(1,k),k=1,na)
325   format('  lower boundary     ',5f10.6)

c   At the positive alternative where power 1-beta is required

      theta=delta
      call gst2(na,inf,theta,zbdy,r,pu,pe,pl,pbya,einf,ierr)
      eratio=einf/ifix
      write(6,331) ierr,pu,pe,pl,pu,einf,eratio
331   format( /'theta=delta:'
     +       //'ierr=',i1
     +       //'Pr{Exit upper boundary}       = ',f10.5
     +        /'Pr{Continue inside boundary}  = ',f10.5
     +        /'Pr{Exit lower boundary}       = ',f10.5
     +       //'Power at theta=delta          = ',f10.5
     +       //'E{Information on termination} = ',f10.5
     +        /'E{Inf}/Fixed sample inf level = ',f10.5)
      write(6,332)
332   format(/'Probability of stopping by analysis'/)
      write(6,333) (pbya(3,k),k=1,na)
333   format('  upper boundary     ',5f10.6)
      write(6,334) (pbya(2,k),k=1,na)
334   format('  between boundaries ',5f10.6)
      write(6,335) (pbya(1,k),k=1,na)
335   format('  lower boundary     ',5f10.6)

      stop
      end
