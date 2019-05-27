
//Code written by Farbod Hassani (12 January 2017)
// Leap-frog -Runge Kutta and Euler method to solve Newtonian Kepler
//Update:16May2019
#include "lib.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "iostream"
#include <fstream>
#include "time.h"
#include <complex>
#include "cmath"
#include "cstdlib"
#include <string>
using namespace std;

int main ()
{
double dt=5.e-2;

  for(int i=1;i<=15;i++)
{

    FILE* Result;
    string base(".txt");
    // string filename;
    std::string filename;
    filename = "./Data_Leap/Newt_leap_data_dt_"+to_string(i)+"_"+base; //s.c_str()
    Result=std::fopen(filename.c_str(),"w");

    //********Initial condition and constants*******//
    double Mearth, M_sun,period, r_per, G_gravity, mercury_eccentricity ,c_light, vmax, Rsch;
    mercury_eccentricity = 0.2056;
    Mearth = 5.972* pow(10.,24.) ;//(*kg*)
    M_sun =1.99*pow(10.,30.)/Mearth ;//(*kg*); (*Normalized to earth's mass*)
    period =87.969*86400/pow(10.,6.) ;//(*Mega second The period of mercury =7.59974 *)
    r_per = 46.001272; //perihelion distance(*m*);//(*R=R(1/10^9)(*(Gm/m)*)(*Gm*);*)
    vmax = 58.98; //
    G_gravity = 6.67*pow(10.,-26.)*Mearth ; //(*(Giga meter)^3/(Mearth (Mega \ second)^2)*);(*G=G(1/10^9)^3(*(Gm/m)^3*)(10^5/1)^2(*(s/d)^2*)(*Gm^3/(kg \ d^2)*)(*d=10^5s*);*)
    c_light = 2.998*pow(10.,8.)*pow(10.,-9.)/pow(10.,-6.);
    vmax = 58.98;//(*(Giga meter)/(Mega second)*)(*vmax=vmax(1/10^9)(*(Gm/m)*)(10^5/1)*)(*(s/d)*)(*Gm/d*)
    Rsch = 2.*G_gravity*M_sun/pow(c_light,2.);
    double cons,x,y,vx,vy,t,tcentury, beta, eccentricity, t_ini, t_fin ;
    ///*********
    //**********
    beta = 1;
    eccentricity = 0.2056;
    t_ini=0;
    t_fin=1.1 * period;
    // t_num=8.e2;
    //***********
    x=r_per * beta *  (1. - eccentricity) /(1. - mercury_eccentricity) ;
    y=0.;
    vx=0.;
    vy= (vmax/sqrt(beta) ) *   sqrt( (1. + eccentricity)/(1. - eccentricity))/ sqrt( (1. + mercury_eccentricity)/(1. - mercury_eccentricity) );



    tcentury=3.65*864;

    cons=G_gravity * M_sun;
      t=0.;
    double r, theta, dr ,r2 ;
      dr=0;
      while  (t<=t_fin)                   // A loop on the time
    {
      //  Er=pow(xr,2.0)+pow(vr,2.0);
      // cout<<"t_f: "<< t_fin << "t" <<t <<endl;
      r=sqrt(pow(x,2.0)+pow(y,2.0));
      theta= atan(y/x);
      if ( std::abs(r_per-x) <200*dt)
      {
        fprintf(Result,"\n  %.17g %.17g %.17g %.17g %.17g %.17g %.17g  ",t , x,y,r ,theta, dr, dt);
      }
      leapfrog(cons,x,y,vx,vy,dt);
      r2=sqrt(pow(x,2.0)+pow(y,2.0));
      dr = r2 - r;
      t=t+dt;
    }
    dt=dt/2;
    cout<<"i "<<i<<endl;
  }


    return 0;
}
