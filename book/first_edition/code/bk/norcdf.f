
c=======================================================================
c
      real*8 function norcdf(x)
c
c-----------------------------------------------------------------------
c
c   NORCDF returns the cumulative distribution of a standard normal
c   distribution.
c
c   This function calls on the NAG library routine  S15ACF.
c   An equivalent routine must be substituted if you do not have access
c   to the NAG library.
c
c                                            C. Jennison, 25 August 1999
c
c-----------------------------------------------------------------------
c
c   INPUT     x                                     REAL*8
c
c   OUTPUT    norcdf    Pr{Z < x} when Z ~ N(0,1)   REAL*8
c
c-----------------------------------------------------------------------
c

      integer ifail
      real*8 x,s15acf

      ifail=0
      norcdf=1-s15acf(x,ifail)

      return
      end
