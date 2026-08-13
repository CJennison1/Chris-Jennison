
c=======================================================================
c
      real*8 function invnor(p)
c
c-----------------------------------------------------------------------
c
c   INVNOR returns a quantile of the standard normal distribution.
c
c   This function calls on the NAG library routine  G01CEF.
c   An equivalent routine must be substituted if you do not have access
c   to the NAG library.
c
c                                            C. Jennison, 25 August 1999
c
c-----------------------------------------------------------------------
c
c   INPUT     p                                        REAL*8
c
c   OUTPUT    invnor    The value of x such that       REAL*8
c                       Pr{Z<x}=p when Z ~ N(0,1)
c
c-----------------------------------------------------------------------
c

      integer ifail
      real*8 p,g01cef

      ifail=0
      invnor=g01cef(p,ifail)

      return
      end
