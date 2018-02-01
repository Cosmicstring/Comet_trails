#ifndef _Generate_SOHO_Image_H_
#define _Generate_SOHO_Image_H_

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define Me 5.97219e24 //kg
#define gama 4.98199486464e-10  //km^3/(kg*d^2)
#define Msol 1.98855e30
#define GMs (gama*Msol)    /**[km³d²]*/
#define AU 1.49597870691e8 // km
//#define n (sqrt(gama*(Msol + Me)/(AU*AU*AU)))
#define n (2*pi/365.25636)
#define T_epoch 2451545.0 // 12:00 1. Januar 2000
#define T_epoch2 2440588.0 // 12:00 1. Januar 1970
#define T_Rosetta_perihel 2457247.58694 
#define T_ISON_perihel 2456625.27616
#define T_Lovejoy_perihel 2455911.51181
#define pi 3.1415926535897932384626433832795
#define DEBUGzaKonverzije 0
#define DEBUGzaPoziciju 0
#define DEBUGzaMisinPredlog 0


typedef struct{
    double x,y,z;
}vectorr;

typedef struct {
  double x, y, z, xd, yd, zd;
} PhaseState;

typedef struct { 
  double a, e, i, longnode, TA, argperi, meananom, meanlongitude;
  double longperi;
} OrbitalElements;

//Function for vectors
vectorr R(vectorr a, vectorr b);
vectorr Norm(vectorr a);
vectorr vp(vectorr a, vectorr b);
vectorr PresekRavniIPrave(vectorr p0, vectorr l, vectorr l0, vectorr N);
double sp(vectorr a, vectorr b);
double intenzitet(vectorr a);

//Some random useful functions
double ConvertToRadians(double x);
double ConvertToDegrees(double x);

//Conversion Functions
double ecc_ano(double e,double l);
double ecc_anohyp(double e,double l);
void keplerian(double GM, PhaseState state,   OrbitalElements * orbel);
void cartesian(double GM, OrbitalElements orbel, PhaseState *state);


#endif