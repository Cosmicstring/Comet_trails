#ifndef _Particle_Orbit_H_
#define _Particle_Orbit_H_

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define DEBUG 1
#define Video 0

const double msol=1.939e30;
const double AU = 1.49597870691e8;
const double gama = 4.98199486464e-10;//6.67e-11 * (86400 * 86400)/(1e9) [km^3/(kg*d^2)]
const double c = (299792.458 * 86400.0); // [km/d]
const double pi = 3.14159265359;
const double L = ((3.839 * (86.4*86.4)*1e21)/2.99792458); // L = L0/c [(kg*km)/d^2]
const double Rs = 6.955e5; //km
const double sat = 0.041666667; // dan
const double minut = 0.000694444; // dan

typedef struct {
 	double x; double y; double z;
} vectorr;

typedef struct{
    double Mass,ro,R;
}meteoroid;

typedef struct {
  double x, y, z, xd, yd, zd;
} PhaseState;

typedef struct {
  double a, e, i, longnode, TA, argperi, meananom;
  vectorr sam;
} OrbitalElements;


// FUNCTIONS WITH VECTORS


inline vectorr Addvectorr(vectorr &a, vectorr &b)
{
    vectorr x;

    x.x = a.x + b.x;
    x.y = a.y + b.y;
    x.z = a.z + b.z;

    return x;
}

inline vectorr vp(vectorr &a, vectorr &b){ // a = r, b = v

vectorr t;
	t.x = b.z*a.y - b.y * a.z;
	t.y = b.x*a.z - b.z * a.x;
	t.z = b.y*a.x - b.x * a.y;
return t;
}

inline double intenzitet(vectorr &a){
	return sqrt(a.x*a.x + a.y*a.y + a.z * a.z);
}

inline double sp(vectorr &a, vectorr &b){

	return (a.x*b.x + a.y*b.y + a.z*b.z);
}

inline double angle(vectorr &a, vectorr &b){

    return sp(a,b)/(intenzitet(a)*intenzitet(b));
}

inline double cosangle(vectorr &a, vectorr &b){

    return sp(a,b)/(intenzitet(a)*intenzitet(b));
}

inline double sinangle(vectorr &a, vectorr &b){

    vectorr t;
    t = vp(a,b);
    return (intenzitet(t)/intenzitet(a)/intenzitet(b));
}

inline vectorr ubrzanja(vectorr &pozicija, vectorr &brzina, double R, double m2){ // Nalazim teta = v/r * t(korak), ort teta = vp(L(ugaoni moment),r(radijus vektor))/ L*r

    vectorr a,teta,Momenti,Vrad;

    double r = intenzitet(pozicija);
    double v = intenzitet(brzina);
    //double beta = 3*L/(16*pi*gama*msol*ro*R);
    double Spr = L/(4*pi*r*r); //L = L0/c
    double S = L/(4*pi);
    double cosalfa = cosangle(pozicija, brzina);
    double sinalfa = sinangle(pozicija,brzina);
    double dteta = v*sinalfa;
    double dr = v*cosalfa;

    Vrad.x = dr*pozicija.x/r;
    Vrad.y = dr*pozicija.y/r;
    Vrad.z = dr*pozicija.z/r;

    double A = pi*R*R;
    double K1 = gama*msol - ((S*A)/m2)*(1 - 2*dr/c);// -> km^3/d^2
    double K3 = ((S*A)/m2)*(1 - 2*dr/c);
    double K2 = (((Spr*A)/c)*dteta)/m2;

    //Check the fraction of radiative and relativistic effects

    /*double RelativisticEffectXclan = (4*gama*msol*pozicija.x/r - (v*cosalfa*v*cosalfa)*pozicija.x + 4.0*sp(pozicija, Vrad)*(Vrad.x * cosalfa));
    double RelativisticEffectYclan = (4*gama*msol*pozicija.y/r - (v*cosalfa*v*cosalfa)*pozicija.y + 4.0*sp(pozicija, Vrad)*(Vrad.y * cosalfa));
    double RelativisticEffectZclan = (4*gama*msol*pozicija.z/r - (v*cosalfa*v*cosalfa)*pozicija.z + 4.0*sp(pozicija, Vrad)*(Vrad.z * cosalfa));
    double RelativisticEffect = (gama*msol/(c*c*r*r*r))*sqrt(RelativisticEffectXclan*RelativisticEffectXclan + RelativisticEffectYclan*RelativisticEffectYclan + RelativisticEffectZclan*RelativisticEffectZclan);
    double GravitacioniClanX = (gama*msol/(r*r*r))*pozicija.x;
    double GravitacioniClanY = (gama*msol/(r*r*r))*pozicija.y;
    double GravitacioniClanZ = (gama*msol/(r*r*r))*pozicija.z;
    double GravitacioniClan = sqrt(GravitacioniClanX*GravitacioniClanX + GravitacioniClanY*GravitacioniClanY + GravitacioniClanZ*GravitacioniClanZ);*/
     /* if( (DEBUG) && ( (Time<=(Tp+0.1)) && (Time>=(Tp-0.1)) )){
        printf("%.10lf\t%.10lf\t%.20lf\n", RelativisticEffect, GravitacioniClan, RelativisticEffect/GravitacioniClan);
        if(MaxEffect<RelativisticEffect/GravitacioniClan) MaxEffect = RelativisticEffect/GravitacioniClan;
    }*/


    double K;

  
    Momenti = vp(pozicija,brzina);
    teta = vp(pozicija,Momenti);

    K = intenzitet(teta);

    a.x = (K1/(r*r*r))*(-pozicija.x) + K2*teta.x/K;
    a.y = (K1/(r*r*r))*(-pozicija.y) + K2*teta.y/K;
    a.z = (K1/(r*r*r))*(-pozicija.z) + K2*teta.z/K;

    return a;
}

/********************************
*********************************/

inline double Max(double r1, double r2)
{
    if(r1 >= r2)
    {
        return r1;
    }
    else return r2;
}

inline double floorfunction(double x,int decp/**na koliko decimala*/)
{
  return floor((x*pow(10,decp)+0.5)/pow(10,decp));
}

inline int trunc_c(double x){

    int t;

    t = (int)x;

    if (t<=x) return t;
    else return(t-1);
}

inline double Adapth(double h, double K, double N, vectorr &r){

    double rp,h1;

    rp = intenzitet(r);

   // printf("rp  %lf\n", rp);

    h1 = sqrt(rp*rp*rp)*K/N;

    //if(h1>0.01) h1 = 0.01;

    return h1;
}

//Conversion Orbit functions

#define PREC_ecc_ano 1e-15  /* no reason that this must be very accurate in code at present */
inline double ecc_ano(double e,double l)
{
    double du,u0,l0;
    du=1.0;
    u0 = l + e*sin(l) + 0.5*e*e*sin(2.0*l);
// also see M+D equation 2.55
                 /* supposed to be good to second order in e, from Brouwer+Clemence
                    u0 is first guess */
    while(fabs(du) > PREC_ecc_ano){
      l0 = u0 - e*sin(u0);
      du = (l - l0)/(1.0 - e*cos(u0));
      u0 += du;  /* this gives a better guess */
// equation 2.58 from M+D
    }
    return u0;
}
// hyperbolic case
inline double ecc_anohyp(double e,double l)
{
    double du,u0,fh,dfh;
    du=1.0;
    u0 = log(2.0*l/e + 1.8); //danby guess
    while(fabs(du) > PREC_ecc_ano){
      fh = e*sinh(u0) -u0 - l;
      dfh = e*cosh(u0) - 1.0;
      du = -fh/dfh;
      u0 += du;
    }
    return u0;
}

void cartesian(double GM, OrbitalElements &orbel, PhaseState *state);
void keplerian(double GM, PhaseState &state,   OrbitalElements *orbel);

inline double Max(double r1, double r2);

#endif