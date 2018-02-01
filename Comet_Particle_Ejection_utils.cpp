#include "Comet_Particle_Ejection.h"
  

  //Mesto za izbacivanje
  ElSpot Spot[MaxN];

/*
*******************************************
****	IMPLEMENTATION OF FUNCTIONS    ****
*******************************************
*/


double intenzitet(vector a)
{
    return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

double ConvertToRadians(double x){
    return(x*pi/180);
}

// Jednacina odrzanja koja se resava (jednacina 3 u radu)
double fun(double x, double r, double sintheta){
    return (A1*x*x*x*x + A2*exp(-6000.0/x)/sqrt(x)-(1361.0/(r*r))*sintheta);
}

double FindTemperature(double r, double sintheta){

    double x1,x2,x,xacc = 1e-6;
    double dx, f, fmid, xmid, rtb;
    int j;

    x1 = 0.01; x2 = 2000;
    f = fun(x1,r,sintheta);
    fmid = fun(x2,r,sintheta);

    x=0.0000052923;
 
    rtb = f < 0.0 ? (dx = x2-x1, x1) : (dx = x1 - x2, x2);
 
    for(j=1; j<100; j++){
        xmid = rtb + (dx *= 0.5);
        fmid = fun(xmid, r, sintheta);
        if(fmid <= 0.0) rtb = xmid;
        if(fabs(dx)<xacc || fmid == 0.0) return rtb;
    }

}

#define PREC_ecc_ano 1e-14  /* no reason that this must be very accurate in code at present */
double ecc_ano(double e,double l)
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
double ecc_anohyp(double e,double l)
{
    double du,u0,fh,dfh;
    du=1.0;
   /* if(l<0){
        while(l<0) l+=2*pi;
    }*/
    //if(DEBUG) printf("DEBUG: %.10lf\t%.10lf\n", l, pow(6*l, 1/3.));
    if((l>=0) && (l<0.1)) u0 = l + (pow(6*l, 1/3.)-l)*e*e;
    else u0 = l + 0.85*e; //danby guess http://research.ijcaonline.org/volume89/number7/pxc3894394.pdf
    while(fabs(du) > PREC_ecc_ano){
      fh = e*sinh(u0) -u0 - l;
      dfh = e*cosh(u0) - 1.0;
      du = -fh/dfh;
      u0 += du;
    }
    return u0;
}


void keplerian(double GM, PhaseState state,   OrbitalElements *orbel)
{
  double rxv_x, rxv_y, rxv_z, hs, h;
  double r, vs, rdotv, rdot, ecostrueanom, esintrueanom, cosnode, sinnode;
  double rcosu, rsinu, u, trueanom, eccanom;

  /* find direction of angular momentum vector */
  rxv_x = state.y * state.zd - state.z * state.yd;
  rxv_y = state.z * state.xd - state.x * state.zd;
  rxv_z = state.x * state.yd - state.y * state.xd;
  hs = rxv_x * rxv_x + rxv_y * rxv_y + rxv_z * rxv_z;
  h = sqrt(hs);

  r = sqrt(state.x * state.x + state.y * state.y + state.z * state.z);
  vs = state.xd * state.xd + state.yd * state.yd + state.zd * state.zd;

  rdotv = state.x * state.xd + state.y * state.yd + state.z * state.zd;
  rdot = rdotv / r;

  orbel->i = acos(rxv_z / h);

  if(rxv_x!=0.0 || rxv_y!=0.0) {
    orbel->longnode = atan2(rxv_x, -rxv_y);
  } else orbel->longnode = 0.0;

  orbel->a = 1.0 / (2.0/r - vs/GM); // could be negative

  ecostrueanom = hs/(GM*r) - 1.0;
  esintrueanom = rdot * h/GM;
  orbel->e = sqrt(ecostrueanom * ecostrueanom + esintrueanom * esintrueanom); /**na pocetku se dobijalo ok*/

  if(esintrueanom!=0.0 || ecostrueanom!=0.0) {
    trueanom = atan2(esintrueanom, ecostrueanom);
  } else trueanom = 0.0;

   orbel->trueanom =trueanom;
  //printf("trueanom: %.10lf\n", trueanom);

  cosnode = cos(orbel->longnode);
  sinnode = sin(orbel->longnode);

  /* u is the argument of latitude */
  rcosu = state.x * cosnode + state.y * sinnode;
  rsinu = (state.y * cosnode - state.x * sinnode)/cos(orbel->i);

  if(rsinu!=0.0 || rcosu!=0.0) {
    u = atan2(rsinu, rcosu);
  } else u = 0.0;

  orbel->argperi = u - trueanom;
  //printf("argperi: %.10lf\n", orbel->argperi);

  double foo = sqrt(fabs(1.0 - orbel->e)/(1.0 + orbel->e));
  if (orbel->e <1.0){
     eccanom = 2.0 * atan(foo*tan(trueanom/2.0));
     orbel->meananom = eccanom - orbel->e * sin(eccanom);
     if (orbel->meananom> M_PI) orbel->meananom-= 2.0*M_PI;
     if (orbel->meananom< -M_PI) orbel->meananom+= 2.0*M_PI;
     // only shift M if elliptic orbit
  }
  else {
     eccanom = 2.0 * atanh(foo*tan(trueanom/2.0));
     orbel->meananom = orbel->e * sinh(eccanom) - eccanom;
  }
 // printf("meananom: %.10lf\n", orbel->meananom);
// printf("M_PI: %.10lf\n", M_PI);
  if (orbel->argperi > M_PI){
    orbel->argperi-= 2.0*M_PI;
   // printf("Upao u prvo if > M_PI: %.10lf\n", orbel->argperi);
  }
  if (orbel->argperi < -M_PI){
    orbel->argperi+= 2.0*M_PI;
    //printf("Upao u drugo if < M_PI %.10lf\n", orbel->argperi);
  }
 /** orbel->i = pi - orbel->i; /** Kada se oduzme od pi, dobije se valjano, ali mora na kraju f-je, jer utice na ostatak parametara
  orbel->longnode = orbel->longnode + pi; /** Izgleda da je potrebno dodati pi nakon njegovog racunanj longnode-a
  orbel->argperi = pi - orbel->argperi;*/
}

void cartesian(double GM, OrbitalElements orbel, PhaseState *state)
{
  double meanmotion, cosE, sinE, foo;
  double x, y, z, xd, yd, zd;
  double xp, yp, zp, xdp, ydp, zdp;
  double cosw, sinw, cosi, sini, cosnode, sinnode;
  double E0,rovera;
  double a = orbel.a;
  double e = orbel.e;
  double i = orbel.i;
  double longnode = orbel.longnode;
  double argperi = orbel.argperi;
  double meananom = orbel.meananom;
  /* double E1, E2, den; */

  /* compute eccentric anomaly */
  if (e<1)
    E0 = ecc_ano(e,meananom);
  else
    E0 = ecc_anohyp(e,meananom);
  // E0 = kepler(e,meananom); // also works

  if (e<1.0){
    cosE = cos(E0);
    sinE = sin(E0);
  }
  else {
    cosE = cosh(E0);
    sinE = sinh(E0);
  }
  a = fabs(a);
  meanmotion = sqrt(GM/(a*a*a));
  foo = sqrt(fabs(1.0 - e*e));
  /* compute unrotated positions and velocities */
  rovera = (1.0 - e * cosE);
  if (e>1.0) rovera *= -1.0;
  x = a * (cosE - e);
  y = foo * a * sinE;

  z = 0.0;
  xd = -a * meanmotion * sinE / rovera;
  yd = foo * a * meanmotion * cosE / rovera;
  zd = 0.0;
  if (e>1.0) x *= -1.0;

  /* rotate by argument of perihelion in orbit plane*/
  cosw = cos(argperi);
  sinw = sin(argperi);

  xp = x * cosw - y * sinw;
  yp = x * sinw + y * cosw;
  zp = z;
  xdp = xd * cosw - yd * sinw;
  ydp = xd * sinw + yd * cosw;
  zdp = zd;

  /* rotate by inclination about x axis */
  cosi = cos(i);
  sini = sin(i);
  x = xp;
  y = yp * cosi - zp * sini;
  z = yp * sini + zp * cosi;
  xd = xdp;
  yd = ydp * cosi - zdp * sini;
  zd = ydp * sini + zdp * cosi;

  /* rotate by longitude of node about z axis */
  cosnode = cos(longnode);
  sinnode = sin(longnode);
  state->x = x * cosnode - y * sinnode;
  state->y = x * sinnode + y * cosnode;
  state->z = z;
  state->xd = xd * cosnode - yd * sinnode;
  state->yd = xd * sinnode + yd * cosnode;
  state->zd = zd;
}

double Adapth(double h, double K, double N, vector r){

    double rp,h1;

    rp = intenzitet(r);
   // printf("rp  %lf\n", rp);

    h1 = sqrt(rp*rp*rp)*K/N;

    return h1;
}

vector AddVector(vector a, vector b){
    vector x;

    x.x = a.x + b.x;
    x.y = a.y + b.y;
    x.z = a.z + b.z;

    return x;
}

double drand ( double low, double high )
{
    return ((double)rand() * (high-low) / (double)RAND_MAX + low); /** formula je OK*/
}

double Convert(double x){

    return ((x*180)/pi);
}

double sp(vector a, vector b){
	return (a.x*b.x + a.y*b.y + a.z*b.z);
}

vector vp(vector a, vector b)
{
    vector t;

	t.x = b.z*a.y - b.y * a.z;
	t.y = b.x*a.z - b.z * a.x;
	t.z = b.y*a.x - b.x * a.y;

    return t;
}

double floorfunction(double x,int decp/**na koliko decimala*/)
{
  return floor((x*pow(10,decp)+0.5)/pow(10,decp));
}

double CalcTemperature(double r, double cosTheta)
{
    return (390*pow(cosTheta,0.25)/sqrt(r));
}

inline int trunc_c(double x){

    int t;

    t = (int)x;

    if (t<=x) return t;
    else return(t-1);
}

double RoundAsSaid(double x)
{
    int t;

    t = trunc_c(x);

    if(fabs(x-t)<=fabs(t-x+1));
}

double Min(double a, double b){
    if(a<b) return a;
    else return b;
}

double Max(double a, double b){

    if(a<b) return b;
    else return a;

}

void GiveMeParticles(FILE*BrzinaiPolozaj, FILE*Vecheck, FILE* Provera, int &MeteorNum, double Zuu, int ni, double Mass, double Rc, double Time, OrbitalElements P, vector r, vector v){

    PhaseState p,P0, Pp;
    OrbitalElements o1, O0;
    int n=0;
    double theta, phi;
    double Ms = 1.989e30; /** [kg]*/
    double s; /**poluprecnik cestice*/
    double Pv,Pv0; /** iz rada 7. **/
    double Z,Z0, Zu; /** formula za izracunavanje iz rada 7.**/
    double mi = (18.02*1.661*1e-27); /** [kg] Whipple explained **/ double mv = 18.02*1e-3; /**kg/mol*/ double mu = 1.660538921e-27; /** u kg **/
    double k = 1.3806488e-23; /** [m2 kg s-2 K-1] **/
    double T;
    double ug;
    double R = 8.3144621;
    double scgs;
    double eta;
    
    double roc = 0.5; /** [kg/m^3] */
    double rocgs = (roc/1000);
    double rop = 2600.0; /** [kg/m^3] */
    double ropcgs =rop/1000.0;
    
    double rp = (intenzitet(r)*1000); /** izrazeno u metrima */
    double rs = (6*AU); /**[m]*/
    double rau = (intenzitet(r)/(AU*1e-3)); /** jer je AU u konstantama gore u metrima */
    
    double C2 = (1.0/(rp*rp) - 1.0/(rs*rs));  double C1,C3;

    double alfaMc, Mcu = 0, Mc, H = 2.838e6; /** [m^2/s^2] */
    double L0 = 3.839e26; /** izrazeno u W */
    double Rkm;
    double rm = intenzitet(r)*1e3;
    double u1,u2,V,x;
    
    double sinTheta, cosTheta, sinPhi, cosPhi, Phi;
    double beta = 0.5; // Negde oko mikrometra
    double alfa = 0.15; /** 15% povrsine komete bi trebalo da bude pokriveno elementarnim celijama, uporediti ovo sa Williamsom*/
    double X = 20;
    double Z00,PV00; /** Produkcija na 1AU */
    double ThetaC, PhiC, SKosinusni, dU, dU00; /** dU je ugaona duzina stranice sfernog kvadrata koji mi je zapravo kosinusni segment */
    
    double MaxReaction = -1;
    
    //Provera za Gravitacione efekte i radijativne efekte uzete formule iz Yeomansa

    double GravitacioniClanX = (GMs*Mass/(rm*rm*rm))*(r.x*1000);
	  double GravitacioniClanY = (GMs*Mass/(rm*rm*rm))*(r.y*1000);
  	double GravitacioniClanZ = (GMs*Mass/(rm*rm*rm))*(r.z*1000);
	  double GravitacioniClan = sqrt(GravitacioniClanX*GravitacioniClanX + GravitacioniClanY*GravitacioniClanY + GravitacioniClanZ*GravitacioniClanZ);
    double ReactiveForceX, ReactiveForceY, ReactiveForceZ, ReactiveForce;
    
    double ra0, rd0;
    double Vispis, Uglic;
    double AUkm = AU*1e-3;
    double RispisAU=0;
    double Ve,Zve,Tve,Pve;
    int i,nx;
    vector ve0, ve, ve2, vi, vi2, k1, l;


    //
    Uglic=10;
    T = FindTemperature(1, 1);
    PV00 = 1.2*exp(-6000.0/T)*1e12;
    Z00 = PV00*sqrt(mi/(2*pi*k*T));

    memset(Spot, 0, sizeof(Spot));
    Zu = 0;

    Rkm = Rc;
    Rc = Rc * 1000; /** posto je Rc u km, a Z je u mol. * m^-2 * s^-1 */

    Mc = (Rc*Rc*L0/(4.0*H))*(1.0/(rm*rm) - 1.0/(36*AU*AU)); /** Prema Williamsovoj formuli */
    Mcu += Mc;
    alfaMc = alfa * Mc;

    // Inicijalni uslovi u 
    T = FindTemperature(rau,1);
    Pv0 = 1.2*exp(-6000.0/T)*1e12; /** N*m^-2 **/
    Z0 = Pv0*sqrt(mi/(2*pi*k*T)); /** Ovo Z0 je m*Z, mass flow rate (Delsemme) */


    //ug = sqrt(8.0*k*T/(pi*mi))/pow(rau,0.25);
    ug = sqrt(8.0*k*T/(pi*mi));
    s = 5.7e-5/(beta*ropcgs);/** Proveriti sa Igorom da li je 0.6 mikrometara [cm] Imam klase cestica 10 mikro, 100, 0.1 mm, 1mm, 1cm */
    s = s*1e-2; /** [m] */

    //Provera za Whipple-ovu brzinu

    if(ni==0){
        for(RispisAU=0.01; RispisAU<10; RispisAU+=0.01){
            Tve = FindTemperature(RispisAU,1);
            Pve = 1.2*exp(-6000.0/Tve)*1e12; /** N*m^-2 **/
            //Z0 = Pv0/sqrt(2*pi*mi*mu*k*T);
            scgs = s*1e2; /** posto je s u cgs gore, rad "Radiation Forces upon spherical particles..", Burns et al.*/
            C1 = scgs*ropcgs*pow(RispisAU,2.25); /** OK su jedinice*/
            C3 = 0.013*rocgs*Rkm; /** Nije ubacena promena radijusa komete u program, da li je uopste uvoditi?*/
            x = (sqrt(1.0/C1 - C3)*sqrt(Rkm)*6.56);
            Zve = Pve*sqrt(mi/(2*pi*k*Tve)); /** Ovo Z0 je m*Z, mass flow rate (Delsemme) */
            Ve = sqrt((Zve*ug)/(2*pi*s*rop*Rc)); /** 2000 metara, 2km poluprecnik ISONa*/
            fprintf(Vecheck,"%.10lf\t%.10lf\t%.10lf\n", RispisAU, Ve, x);
        }
        ni=1;
    }

    ReactiveForceX = (9.0*Z0/Z00)*((r.x*1000)/rm);
    ReactiveForceY = (9.0*Z0/Z00)*((r.y*1000)/rm); /** 9.0 koeficijent sa sajta http://www.researchgate.net/post/How_to_determine_the_trajectory_of_Comet_Ison*/
    ReactiveForceZ = (9.0*Z0/Z00)*((r.z*1000)/rm);
    ReactiveForce = sqrt(ReactiveForceX*ReactiveForceX + ReactiveForceY*ReactiveForceY + ReactiveForceZ*ReactiveForceZ);

    if(MaxReaction < ReactiveForce/GravitacioniClan){ /** Dobija se 1e-11 */
        MaxReaction = ReactiveForce/GravitacioniClan;
      //  printf("%.10lf\t%.10lf\t%.10lf\n", 9.0*Z0/Z00, GMs*Mass/(rm*rm), MaxReaction);
    }

    //Ceo broj cestica se racuna
    if ( fabs(trunc_c(X) - X) <= fabs(trunc_c(X) + 1 - X))
    {
        nx = trunc_c(X);
    }
    else
    {
        nx = trunc_c(X)+1;
    }

    //Zadaje se velicina sfernog kvadrata
    PhiC = drand(0,2*pi); /** Po ekvatoru od 0 do 2pi */
    ThetaC = 2*asin(sqrt(drand(0,1))); /** Od pola do pola */
    SKosinusni = alfa * 4*pi*Rc*Rc;
    dU = sqrt(alfa*4*pi/fabs(cos(ThetaC)));

    while(n<nx)/** Broj cestica je proporcionalan broju 100*(r/r)^-0.24, koji da se proveri sa Igorom! */
    {
        // Uniformno se bira centar sfernog kvadrada
        phi = drand(PhiC-dU/2.0,PhiC+dU/2.0);
        Phi = Convert(phi);
        theta = drand(Max(ThetaC - dU/2.0, 0), Min(ThetaC + dU/2.0, pi)); /** Pogledati kako da stavim drand a da se dobije unifrmno i za mali segment */


        if(sin(theta)!=0){

            u2 = drand(0,1);

            //Trazi se temperatura gasa, pritisak, produkcija i termalna brzina gasa
            T = FindTemperature(rau, sin(theta));
            ug = sqrt(8.0*k*T/(pi*mi)); /** [m/s] */
            Pv = 1.2*exp(-6000.0/T)*1e12;
            Z = Pv*sqrt(mi/(2*pi*k*T)); /** mass flow rate (Delsemme) */

            if(u2<Z/Z0){

                Spot[n].Z = Z;
                Spot[n].phi = phi;
                Spot[n].theta = theta;
                Spot[n].ug = ug;
                Spot[n].s = s;
                Spot[n].ro = rop;

                Zu += Z*SKosinusni/nx;
                Zuu += Zu;

                n++;
            }

            if(n==nx)
            {

                for(i=0; i<n; i++){
                    
                    V = sqrt((Zu*Spot[i].Z/Z0)*Spot[i].ug/(2*pi*Spot[i].s*Spot[i].ro*Rc)); /** 2000 metara, 2km poluprecnik ISONa*/
                    fprintf(Provera, "%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", V, Zu, Spot[i].Z, Spot[i].ug, (Zu*Spot[i].Z/Z0)*Spot[i].ug, 2*pi*s*rop*2000);
          
                    //Zadaje se koordinatni sistem odredjen vektorom momenta impulsa, radijus vektora i tangentne brzine
                    l=vp(r,v);
                    k1=vp(r,l);

                    sinPhi = sin(Spot[i].phi); cosPhi = cos(Spot[i].phi);
                    sinTheta = sin(Spot[i].theta); cosTheta = cos(Spot[i].theta);

                    ve.x = V * ((r.x/intenzitet(r)*sinTheta) + (l.x/intenzitet(l))*cosTheta*sinPhi + k1.x/intenzitet(k1)*cosTheta*cosPhi); /**[m/s]*/
                    ve.y = V * ((r.y/intenzitet(r)*sinTheta) + (l.y/intenzitet(l))*cosTheta*sinPhi + k1.y/intenzitet(k1)*cosTheta*cosPhi); /**[m/s]*/
                    ve.z = V * ((r.z/intenzitet(r)*sinTheta) + (l.z/intenzitet(l))*cosTheta*sinPhi + k1.z/intenzitet(k1)*cosTheta*cosPhi); /**[m/s]*/

                    //Whipple poredjenje uzimam 1 km/s za intenzitet brzine na svim rastojanjima
                    ve2.x = 1000 * ((r.x/intenzitet(r)*sinTheta) + (l.x/intenzitet(l))*cosTheta*sinPhi + k1.x/intenzitet(k1)*cosTheta*cosPhi); /**[m/s]*/
                    ve2.y = 1000 * ((r.y/intenzitet(r)*sinTheta) + (l.y/intenzitet(l))*cosTheta*sinPhi + k1.y/intenzitet(k1)*cosTheta*cosPhi); /**[m/s]*/
                    ve2.z = 1000 * ((r.z/intenzitet(r)*sinTheta) + (l.z/intenzitet(l))*cosTheta*sinPhi + k1.z/intenzitet(k1)*cosTheta*cosPhi); /**[m/s]*/

                    ve0 = ve;

                    double AUkm = AU*1e-3;
                    vector vAU;
                    /*Moraju da budu u AU/dan za Mercury*/

                    ve.x = ve.x*86.4/AUkm;
                    ve.y = ve.y*86.4/AUkm;  /**ok [AU/d]*/
                    ve.z = ve.z*86.4/AUkm;

                    vAU.x = v.x/AUkm;
                    vAU.y = v.y/AUkm;   /**  [AU/d] */
                    vAU.z = v.z/AUkm;

                    vi = AddVector(vAU, ve);

                    ve2.x = ve.x*86.4;
                    ve2.y = ve.y*86.4;  /**ok [km/d]*/
                    ve2.z = ve.z*86.4;

                    vi2 = AddVector(v,ve2);

                    P0.xd = vi.x; P0.yd = vi.y; P0.zd = vi.z;
                    P0.x = r.x; P0.y = r.y; P0.z = r.z;

                    fprintf(BrzinaiPolozaj, " DP%d\t Ep=%.10lf\t d=%.10lf\n %.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%d\t%d\t%d\n", MeteorNum, BeginJD + Time, beta, r.x/AUkm, r.y/AUkm, r.z/AUkm, vi.x, vi.y, vi.z, 0, 0, 0);
                    MeteorNum++;
               }
            }
        }
    }
}
