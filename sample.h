! includes
implicit none
real*16 mu
parameter(mu=4.0q00)
real*16 beta
parameter(beta=1.0q0) !rescales the mercury problem
real*16 dx
parameter(dx=0.0005q0) !grid size
real*16 period
parameter(period=7.59974) !Mercury period in megasecond = 87.96days
real*16 Norbit
parameter(Norbit=3.2) !Number of perihelion we are interested
real*16 tolerance
parameter(tolerance=8.0q-11) !Tolerance  was -11
