c-----------------------------------------------------------------------
c
c   MAINES2  Main program for calculating properties of selected
c   two-sided, error spending tests.
c
c                                          C. Jennison, 9 September 1999
c
c-----------------------------------------------------------------------
c
c   The examples concern tests with error spending functions of the form
c
c     f(t) = min(alpha t^rho,alpha)
c
c   where t = inf/imax
c
c   Constants defining these tests can be found in Table 7.1
c
c   in the book Group Sequential Methods with Applications to Clinical
c   Trials by Jennison & Turnbull.
c
c   The expected sample sizes of the test run under design conditions
c   agree with the values displayed in Table 7.3 (any small difference
c   is because the constant rld has been rounded to 3 decimal places).
c
c   Attained powers when information levels depart from their design
c   values agree with Table 7.4 for the cases "r=1 and pi=0.9 & 1.1"
c   in the notation of that table.
c
c-----------------------------------------------------------------------

      integer na,r,ierr,k
      real*8 inf(100),theta,zbdy(2,100),pu,pe,pl,pbya(3,100),einf,
     +       alpha,beta,delta,rho,ifix,rld,imax,eratio,t1err,
     +       cume(2,100),
     +       invnor

c   Set number of analyses, Type I and II error probabilities, and value
c   of theta at which power 1-beta is required.

      na=5
      alpha=0.05
      beta=0.1
      delta=2.0

      rho=2.0

c   Set the parameter controlling density of grid points in numerical
c   integrations.

      r=16

c   Set constants for the tests.

c   rld is the factor which multiplies the fixed information level to
c   give the target maximum information level.

      rld=1.058

c   Find information level needed by a fixed sample test

      ifix=(invnor(1-alpha/2)+invnor(1-beta))**2/delta**2
      imax=ifix*rld

      write(6,1) alpha,beta,delta,ifix
1     format(/'A fixed sample, one-sided test with Type I'
     +        ' error probability =',f6.3,
     +       /'and Type II error probability =',f6.3,' at theta =',f6.3,
     +        ' requires an'
     +       /'information level =',f6.3,'.')

c   a) Lan & DeMets error spending test when information levels are exactly
c      equal to their design values

      do 110 k=1,na
         inf(k)=(k*imax)/na
         cume(2,k)=0.5*dmin1(alpha*(inf(k)/imax)**rho,alpha)
         cume(1,k)=0.5*dmin1(alpha*(inf(k)/imax)**rho,alpha)
110     continue

      call ersp2(na,inf,cume,r,zbdy,ierr)

c   Under the null hypothesis

      theta=0.0
      call gst2(na,inf,theta,zbdy,r,pu,pe,pl,pbya,einf,ierr)
      eratio=einf/ifix
      t1err=pu+pl
      write(6,121) ierr,pu,pe,pl,t1err,einf,eratio
121   format(//'Two-sided error spending test'
     +        /'-----------------------------'
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

c   b) Maximum duration error spending test when information levels are
c      a little smaller than their design values.
c
c      Note that the information sequence does not reach imax but, as
c      specified in the definition of a "maximum duration" test, the
c      Type I error is brought up to its overall value at the final
c      analysis.

      do 210 k=1,na
         inf(k)=0.9*(k*imax)/na
         cume(2,k)=0.5*dmin1(alpha*(inf(k)/imax)**rho,alpha)
         cume(1,k)=0.5*dmin1(alpha*(inf(k)/imax)**rho,alpha)
         if(k.eq.na) then
            cume(2,k)=0.5*alpha
            cume(1,k)=0.5*alpha
         endif
210     continue

      call ersp2(na,inf,cume,r,zbdy,ierr)

c   Under the null hypothesis

      theta=0.0
      call gst2(na,inf,theta,zbdy,r,pu,pe,pl,pbya,einf,ierr)
      eratio=einf/ifix
      t1err=pu+pl
      write(6,221) ierr,pu,pe,pl,t1err,einf,eratio
221   format(//'Two-sided error spending test, low inf. levels'
     +        /'----------------------------------------------'
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

c   c) Maximum duration or maximum information error spending test when
c      information levels are a little larger than their design values.
c
c      Since the information level imax is reached, the same test results
c      in either a "maximum duration" or "maximum information" formulation.

      do 310 k=1,na
         inf(k)=1.1*(k*imax)/na
         cume(2,k)=0.5*dmin1(alpha*(inf(k)/imax)**rho,alpha)
         cume(1,k)=0.5*dmin1(alpha*(inf(k)/imax)**rho,alpha)
310     continue

      call ersp2(na,inf,cume,r,zbdy,ierr)

c   Under the null hypothesis

      theta=0.0
      call gst2(na,inf,theta,zbdy,r,pu,pe,pl,pbya,einf,ierr)
      eratio=einf/ifix
      t1err=pu+pl
      write(6,321) ierr,pu,pe,pl,t1err,einf,eratio
321   format(//'Two-sided error spending test, high inf. levels'
     +        /'-----------------------------------------------'
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
