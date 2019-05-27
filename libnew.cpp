 // *This subroutine is written by Farbod Hassani (Jan 2017)
//*Last update :18 January 2017 (Farbod Hassani)



#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "iostream"
#include "stdlib.h"
#include "time.h"
#include <complex>
#include "cmath"
#include "cstdlib"

using namespace std;






//------------ Constant Definition ----------------

#define PI 3.1415926535897932384626433832795028841971693993751

//------------Sub-Routines/ Functions-------------------

//*************************************


// ***************MATHEMATHICAL FUNCTIONS *******************

//--------Sub-Routine for squared ----------------

double square (double x )
{
    double y;
    y=x*x;

    return (y);
}


//---------------------------------------------------


double f_funcy(double cons,double x,double y){
    double rr,coef;
    coef=pow(square(x)+square(y),1.5);
    rr=-cons*y/coef;
    return (rr);
}
double f_funcx(double cons,double x,double y){
    double rr,coef;
    coef=pow(square(x)+square(y),1.5);
    rr=-cons*x/coef;
    return (rr);
}
double v_func(double v){
    return (v);
}

void leapfrog (double& cons,double& x,double& y,double& vx,double& vy, double& h){
    // kick-drift-kick
    vx = vx + h * f_funcx(cons,x,y) /2.;
    vy = vy + h * f_funcy(cons,x,y) /2.;
    x =  x  + h * vx ;
    y =  y  + h * vy;
    vx = vx + h * f_funcx(cons,x,y) /2.;
    vy = vy + h * f_funcy(cons,x,y) /2.;
}

void Euler (double& cons,double& x,double& y, double& vx, double& vy, double& h)
{
    vx = vx + h * f_funcx(cons,x,y);
    vy = vy + h * f_funcy(cons,x,y) ;
    x =  x  + h * vx ;
    y =  y  + h * vy;
}

void Mid_point_RK2 (double& cons,double& x,double& y,double& vx,double& vy, double& h){
    double kx1,ky1,kx2,ky2,kvx1,kvx2,kvy1,kvy2;
    kvx1 = f_funcx(cons,x,y);
    kvy1 = f_funcy(cons,x,y);
    kx1  = v_func(vx);
    ky1  = v_func(vy);

    kvx2 = f_funcx(cons,x+kx1*h/2,y+ky1*h/2);
    kvy2 = f_funcy(cons,x+kx1*h/2,y+ky1*h/2);
    kx2  = v_func(vx+kvx1*h/2);
    ky2  = v_func(vy+kvy1*h/2);


    vx = vx+ h*kvx2;
    vy = vy+ h*kvy2;
    x  = x + h*kx2;
    y  = y + h*ky2;
}


void RK (double& cons,double& x,double& y,double& vx,double& vy, double& h){
    double kx1,ky1,kx2,ky2,kx3,ky3,kx4,ky4,kvx1,kvx2,kvx3,kvx4,kvy1,kvy2,kvy3,kvy4;
    kvx1 = f_funcx(cons,x,y);
    kvy1 = f_funcy(cons,x,y);
    kx1  = v_func(vx);
    ky1  = v_func(vy);

    kvx2 = f_funcx(cons,x+kx1*h/2,y+ky1*h/2);
    kvy2 = f_funcy(cons,x+kx1*h/2,y+ky1*h/2);
    kx2  = v_func(vx+kvx1*h/2);
    ky2  = v_func(vy+kvy1*h/2);


    kvx3 = f_funcx(cons,x+kx2*h/2,y+ky2*h/2);
    kvy3 = f_funcy(cons,x+kx2*h/2,y+ky2*h/2);
    kx3  = v_func(vx+kvx2*h/2);
    ky3  = v_func(vy+kvy2*h/2);

    kvx4 = f_funcx(cons,x+kx3*h,y+ky3*h);
    kvy4 = f_funcy(cons,x+kx3*h,y+ky3*h);
    kx4  = v_func(vx+kvx3*h);
    ky4  = v_func(vy+kvy3*h);

    vx = vx+h*(kvx1+2*kvx2+2*kvx3+kvx4)/6.;
    vy = vy+h*(kvy1+2*kvy2+2*kvy3+kvy4)/6.;
    x  = x+h*(kx1+2*kx2+2*kx3+kx4)/6.;
    y  = y+h*(ky1+2*ky2+2*ky3+ky4)/6.;
}
