c-----------------------------------------------------------------------
c
c   MAIN1  Main program for calculating properties of selected one-sided
c   group sequential tests.
c
c                                            C. Jennison, 25 August 1999
c
c-----------------------------------------------------------------------
c
c   Two examples are presented
c
c     a)  a Pampallona & Tsiatis test with 4 analyses,
c         alpha=0.05, beta=0.2 and shape parameter del=0.25,
c         and power 1-beta to be achieved at theta=delta=1.0.
c
c     a)  a Pampallona & Tsiatis test with 5 analyses,
c         alpha=0.05, beta=0.1 and shape parameter del=-0.5,
c         and power 1-beta to be achieved at theta=delta=0.2.
c
c   Constants defining these tests are taken from Tables 4.1 and 4.2.
c   in the book Group Sequential Methods with Applications to Clinical
c   Trials by Jennison & Turnbull.
c
c   The expected sample sizes of these tests agree with the values
c   displayed in the same tables.  (Note that error probabilities do
c   differ a little from specified values as constants ct1, ct2 and rt
c   have been rounded to 3 decimal places; similarly, E{inf} may differ
c   in the last decimal place from the tabulated value which is based
c   on slightly more accurate values of ct1, ct2 and rt.)
c
c-----------------------------------------------------------------------

      integer na,r,ierr,k
      real*8 inf(100),theta,zbdy(2,100),pu,pl,pbya(2,100),einf,
     +       alpha,beta,delta,ifix,ct1,ct2,rt,del,eratio,
     +       invnor

c   Example (a)
c   ===========

c   Set number of analyses, Type I and II error probabilities, and value
c   of theta at which power 1-beta is required.

      na=4
      alpha=0.05
      beta=0.2
      delta=2.0

c   Set the parameter controlling density of grid points in numerical
c   integrations.

      r=16

c   Set constants for the Pampallona & Tsiatis test.  Here, ct1 and
c   ct2 denote the constants tilde{c}_1 and tilde{c}_2 in Jennison &
c   Turnbull's notation.

      del=0.25

      ct1=1.741
      ct2=1.093
      rt=1.299

c   Find information level needed by a fixed sample one-sided test

      ifix=(invnor(1-alpha)+invnor(1-beta))**2/delta**2

      write(6,101) del,alpha,beta,delta,ifix
101   format(/'Pampallona and Tsiatis test, Del =',f6.2
     +        /'----------------------------------------',
     +       //'A fixed sample, one-sided test with Type I'
     +        ' error probability =',f6.3,
     +       /'and Type II error probability =',f6.3,' at theta =',f6.3,
     +        ' requires an'
     +       /'information level =',f8.3,'.')

      do 110 k=1,na
         inf(k)=(k*rt*ifix)/na

c   Using formula (4.3) of the book would give
c
c        zbdy(2,k)=ct1*(k*1.0d0/na)**(del-0.5)
c        zbdy(1,k)=delta*dsqrt(inf(k))-ct2*(k*1.0d0/na)**(del-0.5)
c
c   However, because of rounding in ct1, ct2 and rt2, this formula
c   does not necessarily give zbdy(1,na)=zbdy(2,na) to very high
c   accuracy.  It is, therefore, better to use formula (4.5), from
c   which:

         zbdy(2,k)=ct1*(k*1.0d0/na)**(del-0.5)
         zbdy(1,k)=(ct1+ct2)*dsqrt(k*1.0d0/na)
     +             -ct2*(k*1.0d0/na)**(del-0.5)

110    continue

c   Under the null hypothesis

      theta=0.0
      call gst1(na,inf,theta,zbdy,r,pu,pl,pbya,einf,ierr)
      eratio=einf/ifix
      write(6,121) ierr,pu,pl,pu,einf,eratio
121   format( /'theta=0:'
     +       //'ierr=',i1
     +       //'Pr{Exit upper boundary}       = ',f10.5
     +        /'Pr{Exit lower boundary}       = ',f10.5
     +       //'Type I error probability      = ',f10.5
     +       //'E{Information on termination} = ',f10.5
     +        /'E{Inf}/Fixed sample inf level = ',f10.5)
      write(6,122)
122   format(/'Probability of stopping by analysis'/)
      write(6,123) (pbya(2,k),k=1,na)
123   format('  upper boundary     ',5f10.6)
      write(6,124) (pbya(1,k),k=1,na)
124   format('  lower boundary     ',5f10.6)

c   At the positive alternative where power 1-beta is required

      theta=delta
      call gst1(na,inf,theta,zbdy,r,pu,pl,pbya,einf,ierr)
      eratio=einf/ifix
      write(6,131) ierr,pu,pl,pu,einf,eratio
131   format( /'theta=delta:'
     +       //'ierr=',i1
     +       //'Pr{Exit upper boundary}       = ',f10.5
     +        /'Pr{Exit lower boundary}       = ',f10.5
     +       //'Power at theta=delta          = ',f10.5
     +       //'E{Information on termination} = ',f10.5
     +        /'E{Inf}/Fixed sample inf level = ',f10.5)
      write(6,132)
132   format(/'Probability of stopping by analysis'/)
      write(6,133) (pbya(2,k),k=1,na)
133   format('  upper boundary     ',5f10.6)
      write(6,134) (pbya(1,k),k=1,na)
134   format('  lower boundary     ',5f10.6)

c   Example (b)
c   ===========

c   Set number of analyses, Type I and II error probabilities, and value
c   of theta at which power 1-beta is required.

      na=5
      alpha=0.05
      beta=0.1
      delta=0.2

c   Set the parameter controlling density of grid points in numerical
c   integrations.

      r=16

c   Set constants for the Pampallona & Tsiatis test.  Here, ct1 and
c   ct2 denote the constants tilde{c}_1 and tilde{c}_2 in Jennison &
c   Turnbull's notation.

      del=-0.5

      ct1=1.648
      ct2=1.320
      rt=1.029

c   Find information level needed by a fixed sample one-sided test

      ifix=(invnor(1-alpha)+invnor(1-beta))**2/delta**2

      write(6,201) del,alpha,beta,delta,ifix
201   format(//'Pampallona and Tsiatis test, Del =',f6.2
     +        /'----------------------------------------',
     +       //'A fixed sample, one-sided test with Type I'
     +        ' error probability =',f6.3,
     +       /'and Type II error probability =',f6.3,' at theta =',f6.3,
     +        ' requires an'
     +       /'information level =',f8.3,'.')

      do 210 k=1,na
         inf(k)=(k*rt*ifix)/na

c   Using formula (4.3) of the book would give
c
c        zbdy(2,k)=ct1*(k*1.0d0/na)**(del-0.5)
c        zbdy(1,k)=delta*dsqrt(inf(k))-ct2*(k*1.0d0/na)**(del-0.5)
c
c   However, because of rounding in ct1, ct2 and rt2, this formula
c   does not necessarily give zbdy(1,na)=zbdy(2,na) to very high
c   accuracy.  It is, therefore, better to use formula (4.5), from
c   which:

         zbdy(2,k)=ct1*(k*1.0d0/na)**(del-0.5)
         zbdy(1,k)=(ct1+ct2)*dsqrt(k*1.0d0/na)
     +             -ct2*(k*1.0d0/na)**(del-0.5)

210    continue

c   Under the null hypothesis

      theta=0.0
      call gst1(na,inf,theta,zbdy,r,pu,pl,pbya,einf,ierr)
      eratio=einf/ifix
      write(6,221) ierr,pu,pl,pu,einf,eratio
221   format( /'theta=0:'
     +       //'ierr=',i1
     +       //'Pr{Exit upper boundary}       = ',f10.5
     +        /'Pr{Exit lower boundary}       = ',f10.5
     +       //'Type I error probability      = ',f10.5
     +       //'E{Information on termination} = ',f10.5
     +        /'E{Inf}/Fixed sample inf level = ',f10.5)
      write(6,222)
222   format(/'Probability of stopping by analysis'/)
      write(6,223) (pbya(2,k),k=1,na)
223   format('  upper boundary     ',5f10.6)
      write(6,224) (pbya(1,k),k=1,na)
224   format('  lower boundary     ',5f10.6)

c   At the positive alternative where power 1-beta is required

      theta=delta
      call gst1(na,inf,theta,zbdy,r,pu,pl,pbya,einf,ierr)
      eratio=einf/ifix
      write(6,231) ierr,pu,pl,pu,einf,eratio
231   format( /'theta=delta:'
     +       //'ierr=',i1
     +       //'Pr{Exit upper boundary}       = ',f10.5
     +        /'Pr{Exit lower boundary}       = ',f10.5
     +       //'Power at theta=delta          = ',f10.5
     +       //'E{Information on termination} = ',f10.5
     +        /'E{Inf}/Fixed sample inf level = ',f10.5)
      write(6,232)
232   format(/'Probability of stopping by analysis'/)
      write(6,233) (pbya(2,k),k=1,na)
233   format('  upper boundary     ',5f10.6)
      write(6,234) (pbya(1,k),k=1,na)
234   format('  lower boundary     ',5f10.6)

      stop
      end
