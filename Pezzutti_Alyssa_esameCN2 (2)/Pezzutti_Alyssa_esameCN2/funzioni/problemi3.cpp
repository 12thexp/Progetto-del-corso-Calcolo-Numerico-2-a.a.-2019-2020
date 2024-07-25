#include <math.h>
#include <iostream>
#include "/home/squiddy/Downloads/Calcolo_2_LAB/funzioni/problemi.h"

using namespace std;
typedef double Real;
const Real pi=4.0*atan(1.0);


Real l1 = 1.0; Real l2 = 1.0;
Real m1 = 1.0; Real m2 = 0.99;
Real r = 0.5;
Real c0 = 0.5;
Real ft = 0;

Real C1 = c0*r*r/(m1*l1*l1);
Real C2 = c0*r*r/(m2*l2*l2);

void f_pendolo(Real F[],Real t,Real U[])    //t tempo iniziale, anche se l'eq non dipende da t lo devo passare
{
///U[0]=φ1, U[1]=φ2, U[2]=φ1', U[3]=φ2'
///F[0]=φ1', F[1]=φ2', F[2]=φ1", F[3]=φ2"
    if (abs(t-1) <= 1)
        ft = sqrt(1-pow(1-t,2));
    else
        ft = 0;

    F[0] = U[2];
    F[1] = U[3];
    F[2] = -(sin(U[0])/l1) - (c0*r*r/(m1*l1*l1))*(sin(U[0])-sin(U[1]))*cos(U[0]) + ft;
    F[3] = -(sin(U[1])/l2) - (c0*r*r/(m2*l2*l2))*(sin(U[1])-sin(U[0]))*cos(U[1]);
}
//// definisco dati iniziali
void dati_pendolo(Real *t0,Real *T,Real v[])
{
    *t0 = 0.0;
    *T = 250.0;
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    v[3] = 0.0;
}

void jac_pendolo(Real Jf[], Real t, Real U[])
{
    Jf[0] = 0.0;
    Jf[1] = 0.0;
    Jf[2] = 1.0;
    Jf[3] = 0.0;

    Jf[4] = -(1/l1)*cos(U[0]) - C1*(pow(cos(U[0]),2) - pow(sin(U[0]),2) + sin(U[1])*sin(U[0]));
    Jf[5] = C1*cos(U[1])*cos(U[0]);
    Jf[6] = 0.0;
    Jf[7] = 0.0;

    Jf[8] = 0.0;
    Jf[9] = 0.0;
    Jf[10] = 0.0;
    Jf[11] = 1.0;

    Jf[12] = C2*cos(U[1])*cos(U[0]);
    Jf[13] = -(1/l2)*cos(U[1]) - C2*(pow(cos(U[1]),2) - pow(sin(U[1]),2) + sin(U[1])*sin(U[0]));
    Jf[14] = 0.0;
    Jf[15] = 0.0;

}

