#include <stdio.h>
#include <math.h>
#define X0 1 // initial coordinate
#define V0 0 // initial velocity
#define T0 0 // initial time
#define TFAC 10 // number of time steps per oscillation period
#define NSTEPS (TFAC * 10) // simulate for ten periods
#define TWOPI (2*acos(-1.0)) //  2 pi
#define MASS 1
#define K 1

void main()
// Integration of the simple harmonic oscillator with
//   the Velocity-Verlet algorithm. 

{

double x, v, t; // coordinate and velocity and time
double Tper; // period of oscillation
double omega; // angular frequency 
int i;
double Etot, Ekin, Epot, Etot_shadow; //energies
double deltat; // time step
double atmp1,atmp2;


// initialize
x = X0;
v = V0;
t = T0;

// compute properties
omega = sqrt(K/MASS);
Tper = TWOPI / omega;
deltat = Tper / TFAC;

// write out initial point and energy
Ekin = 0.5 * MASS * v * v;
Epot = 0.5 * K * x * x;
Etot = Ekin + Epot;
Etot_shadow = Ekin / (1.0 - (omega * deltat / 2.0) * (omega * deltat / 2.0) ) + Epot;
printf("# time/T            x            p          Etot\n");
printf(" %12.6f %12.6f %12.6f %12.6f\n",t/Tper,x,v*MASS,Etot);

/* Velocity-Verlet algorithm for the position and velocity */

for (i = 0; i < NSTEPS; ++i)
  {
// eqs of motion are
// xdot = p/m = velocity
// pdot = -kx = acceleration * mass
// (or vdot = -kx / MASS ) 

// acceleration at this position 
  atmp1 = - K * x / MASS ;
// get new position 
  x +=  deltat * v + 0.5 * deltat * deltat * atmp1 ;
// acceleration at new position 
  atmp2 = - K * x / MASS ;
// get new velocity 
  v += 0.5 * deltat * (atmp1 + atmp2);
 
// energies  
Ekin = 0.5 * MASS * v * v;
Epot = 0.5 * K * x * x;
Etot = Ekin + Epot;
//Etot_shadow = Ekin / (1.0 - (omega * deltat / 2.0) * (omega * deltat / 2.0) ) + Epot;
// increment time 
t += deltat;
printf(" %12.6f %12.6f %12.6f %12.6f\n",t/Tper,x,v*MASS,Etot);
  }
}
