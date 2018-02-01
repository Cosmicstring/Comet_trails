#include "Particle_Orbit.h"

/*
*******************************************
****    IMPLEMENTATION OF FUNCTIONS    ****
*******************************************
*/


void keplerian(double GM, PhaseState &state,   OrbitalElements *orbel)
{
  double rxv_x, rxv_y, rxv_z, hs, h;
  double r, vs, rdotv, rdot, ecostrueanom, esintrueanom, cosnode, sinnode;
  double rcosu, rsinu, u, trueanom, eccanom;

  /* find direction of angular momentum vectorr */
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

  orbel->TA = trueanom;

  cosnode = cos(orbel->longnode);
  sinnode = sin(orbel->longnode);

  /* u is the argument of latitude */
  rcosu = state.x * cosnode + state.y * sinnode;
  rsinu = (state.y * cosnode - state.x * sinnode)/cos(orbel->i);

  if(rsinu!=0.0 || rcosu!=0.0) {
    u = atan2(rsinu, rcosu);
  } else u = 0.0;

  orbel->argperi = u - trueanom;
 
  double foo = sqrt(fabs(1.0 - orbel->e)/(1.0 + orbel->e));
  if (orbel->e <1.0){
     eccanom = 2.0 * atan(foo*tan(trueanom/2.0));
     orbel->meananom = eccanom - orbel->e * sin(eccanom);
     if (orbel->meananom> M_PI) orbel->meananom-= 2.0*M_PI;
     if (orbel->meananom< -M_PI) orbel->meananom+= 2.0*M_PI;
  }
  else {
     eccanom = 2.0 * atanh(foo*tan(trueanom/2.0));
     orbel->meananom = orbel->e * sinh(eccanom) - eccanom;
  }

  if (orbel->argperi > M_PI){
    orbel->argperi-= 2.0*M_PI;
  }
  if (orbel->argperi < -M_PI){
    orbel->argperi+= 2.0*M_PI;
  }
}

void cartesian(double GM, OrbitalElements &orbel, PhaseState *state)
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