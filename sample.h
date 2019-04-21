! includes
implicit none
real*16 mu
parameter(mu=4.0q00)
real*16 beta
parameter(beta=1.0q0) !rescales the mercury problem
real*16 dx
parameter(dx=0.001q0) !grid size
real*16 period
parameter(period=7.59974) !Mercury period in megasecond = 87.96days
real*16 Norbit
parameter(Norbit=0.2) !Number of perihelion we are interested
real*16 tolerance
parameter(tolerance=1.0q-11) !Tolerance
real*16 theta
parameter(theta= 0*3.1415926535897932384626433832795028841971693993751/9)
real*16 eccentricity
parameter(eccentricity=0.55) !For 0 we recover the mercury and sun orbit eccentricity
