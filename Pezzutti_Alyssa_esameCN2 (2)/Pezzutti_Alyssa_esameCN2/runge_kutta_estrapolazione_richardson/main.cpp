#include <iostream>
#include <math.h>
#include <fstream>
#include"../../funzioni/mat_vet.h"
#include "../../funzioni/stampe.h"
#include "../../funzioni/problemi.h"
#include "../../funzioni/Runge_Kutta_matrici.h"
using namespace std;
typedef double Real;


int main () {

    // predispongo file stampa
    char n_file[21]= {0};
    cout << "dammi nome file di stampa(max 20 caratteri)";
    cin >> n_file ;
    ofstream prt(n_file);
    prt.precision(20);
    // definisco puntatori al tipo di funzioni che definiscono il problema
    void(*effe)(Real *,Real,Real *);
    void(*dati)(Real *,Real *,Real *);
    //
    // dati del problema **********************************   //le uniche righe che devo cambiare
    const int d=4;
    effe=f_pendolo;
    dati=dati_pendolo;
    //
    // ******************************************
    // inizializzo
    Real t,T,u[d],U2[d],U1[d], fc;   //in u sol. approssimata, fc = fattore correzione
    dati(&t,&T,u);
    // scelta del metodo *****************************
    char metodo[]="RK4";// scelta del metodo RK
    const int s=4; // numero degli stadi del metodo scelto;
    const Real p=4.0; //ordine del metodo
    int f_cnt = 0;

    // dimensiono array di Butcher *********************************
    Real b[s]={0};
    Real c[s]={0};
    Real A[s][s]={{0}}; //
    // assegno valori array di Butcher
    set_Runge_Kutta_esp(metodo,A[0],b,c);

    // N e tau
    unsigned long N;
    cout << "inserisci numero di passi N  = ";
    cin >> N;
    Real tau=(T-t)/Real(N);

    // inizio il ciclo sul tempo
    stampa_estr(d,t,u,0,&prt);   //stampo valori iniziali
    stampa_estr(d,t,u,0,0);



    for (unsigned long n=1; n<=N; n++) {

        ///calcolo U2 = U(1/2)
        Real K2[s][d] = {{0}};
        effe(K2[0],t,u);
        f_cnt++;
        for (int i=1; i<s; i++) {
            Real V2[d] = {0};
            matmat(A[i],K2[0],V2,1,i,d);
            for (int k=0; k<d; k++) {
                V2[k] = u[k] + (tau/2)*V2[k];
            }
            effe(K2[i],t+c[i]*(tau/2),V2);
            f_cnt++;
        }
        Real aux2[d] = {0};
        matmat(b,K2[0],aux2,1,s,d);

        for (int k=0; k<d; k++) {
            U2[k] = u[k] + (tau/2)*aux2[k];
        }


        ///calcolo U1
        Real K1[s][d] = {{0}};
        effe(K1[0],t+tau/2,U2);
        f_cnt++;
        for (int i=1; i<s; i++) {
            Real V1[d] = {0};
            matmat(A[i],K1[0],V1,1,i,d);
            for (int k=0; k<d; k++) {
                V1[k] = U2[k] + (tau/2)*V1[k];
            }
            effe(K1[i],t+tau/2+c[i]*(tau/2),V1);
            f_cnt++;
        }
        Real aux1[d] = {0};
        matmat(b,K1[0],aux1,1,s,d);

        for (int k=0; k<d; k++) {
            U1[k] = U2[k] + (tau/2)*aux1[k];
        }



        ///calcolo U^, il fattore di correzione e U(n+1)
        Real K[s][d] = {{0}};      //calcolo i Ki iniziali
        effe(K[0],t,u);
        f_cnt++;

        for (int i=1; i<s; i++) {
            Real V[d] = {0};
            matmat(A[i],K[0],V,1,i,d);          // V = aij*Ki
            for (int k=0; k<d; k++) {
                V[k] = u[k] + tau*V[k];             // V = Un + tau*Σ(aij*Kj)
            }
            effe(K[i],t+c[i]*tau,V);           //calcolo i Ki
            f_cnt++;
        }
        Real aux[d] = {0};
        matmat(b,K[0],aux,1,s,d);       //aux = bi*Ki

        Real dU[d] = {0};
        for (int k=0; k<d; k++) {
            u[k] = u[k] + tau*aux[k];  ///U^

            dU[k] = U1[k]-u[k];

            u[k] = U1[k] + (dU[k])/(pow(2.0,p)-1.0);   ///U(n+1)
        }

        Real norm = 0;                  ///fattore di correzione
        for (int j=0; j<d; j++) {
            norm += dU[j]*dU[j];
        }
        fc = sqrt(norm) / (pow(2.0,p)-1.0);


//        fc = norm2(dU,d)/(pow(2.0,p)-1);    ///fattore di correzione

    t+=tau;

    stampa_estr(d,t,u,fc,&prt);
    stampa_estr(d,t,u,fc,0);

    }
    cout << "cnt f = " << f_cnt << endl;
    return 0;
}
