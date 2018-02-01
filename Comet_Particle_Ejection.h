#ifndef _Comet_Particle_Ejection_H_
#define _Comet_Particle_Ejection_H_

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
/**sa sajta http://jenab6.livejournal.com/34132.html*/



#define DEBUG 0

#define AU 1.49597870691e11    /**[m]*/
#define GMs 1.32712440018e20    /**[m³s²]*/
#define pi 3.1415926535897932384626433832795
#define msol 1.939e30 /**[kg]*/
#define gama 6.67384e-11 /**m^3 kg^-1 s^-2*/
#define Rs 6.955e8 /** [m] */
#define A1 5.67037321e-8
#define Na 6.0221412927e23
#define A2 (78.6975555e39/Na)
#define BeginJD 2456882.5
#define MaxN 100

typedef struct
{
    double x,y,z;
}vector;

typedef struct
{
    double Z,ug,s,ro,theta,phi;
}ElSpot;

typedef struct
{
    double Mass,R,ro;
    double a, e, i, longnode, TA, argperi, meananom, trueanom;
    vector sam;
}comet;

typedef struct {
  double x, y, z, xd, yd, zd;
} PhaseState;

typedef struct {
  double a, e, i, longnode, TA, argperi, meananom, trueanom;
  vector sam;
} OrbitalElements;
    

//Functions for vectors
vector vp(vector a, vector b);
vector AddVector(vector a, vector b);
double sp(vector a, vector b);
double intenzitet(vector a);

//Temperature
double CalcTemperature(double r, double cosTheta);
double fun(double x, double r, double sintheta);
double FindTemperature(double r, double sintheta);

//Helpful functions
double Convert(double x);
double drand ( double low, double high );
double floorfunction(double x,int decp/**na koliko decimala*/);
double RoundAsSaid(double x);
double ConvertToRadians(double x);
int trunc_c(double x);

//Integration Functions
double Adapth(double h, double K, double N, vector r);
void GenMat(double **GridT, double r);

//Conversion functions (they take radians)
double ecc_ano(double e,double l);
double ecc_anohyp(double e,double l);
void keplerian(double GM, PhaseState state,   OrbitalElements * orbel);
void cartesian(double GM, OrbitalElements orbel, PhaseState *state);

//Generate particles
void GiveMeParticles(FILE* BrzinaoPolozaj, FILE* Vecheck, FILE* Provera, int &MeteorNum, double Zuu, int ni, double Mass, double Rc, double Time, OrbitalElements p, vector r, vector v);

#endif