      real*16 function zeroin(ax,bx,f,tol)
      implicit real*16 (a-h,o-z)
      external f
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result ( .ge. 0.0)
c
c
c  output..
c
c  zeroin abcissa approximating a zero of  f  in the interval ax,bx
c
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
c  is the relative machine precision.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
c
c
c  compute eps, the relative machine precision
c
      eps = 1.0Q0
   10 eps = eps/2.0Q0
      tol1 = 1.0Q0 + eps
      if (tol1 .gt. 1.0Q0) go to 10
c initialization
c
      a = ax
      b = bx
      fa = f(a)
      fb = f(b)

      !print *,f(a),f(b),'zero'
      if(f(a)*f(b)>0)then
      !   print *,a,b,'no sign change'
         stop
      endif
c
c begin step
c
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if (abs(fc) .ge. abs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c convergence test
c
 40   continue
      tol1 = 2.0Q0*eps*abs(b) + 0.5Q0*tol
      xm = .5Q0*(c - b)
      if (abs(xm) .le. tol1) go to 90
      if (fb .eq. 0.0Q0) go to 90
c
c is bisection necessary
c
      if (abs(e) .lt. tol1) go to 70
      if (abs(fa) .le. abs(fb)) go to 70
c
c is quadratic interpolation possible
c
      if (a .ne. c) go to 50
c
c linear interpolation
c
      s = fb/fa
      p = 2.0Q0*xm*s
      q = 1.0Q0 - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0Q0*xm*q*(q - r) - (b - a)*(r - 1.0Q0))
      q = (q - 1.0Q0)*(r - 1.0Q0)*(s - 1.0Q0)
c
c adjust signs
c
   60 if (p .gt. 0.0Q0) q = -q
      p = abs(p)
c
c is interpolation acceptable
c
      if ((2.0Q0*p) .ge. (3.0Q0*xm*q - abs(tol1*q))) go to 70
      if (p .ge. abs(0.5Q0*e*q)) go to 70
      e = d
      d = p/q
      go to 80
c
c bisection
c
   70 d = xm
      e = d
c
c complete step
c
   80 a = b
      fa = fb
      if (abs(d) .gt. tol1) b = b + d
      if (abs(d) .le. tol1) b = b + sign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/abs(fc))) .gt. 0.0Q0) go to 20
      go to 30
c
c done
c
   90 zeroin = b
      return
      end
