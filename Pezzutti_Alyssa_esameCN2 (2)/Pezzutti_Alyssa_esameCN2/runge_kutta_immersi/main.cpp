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
    char n_file[21] = {0};
    cout << "dammi nome file di stampa(max 20 caratteri)";
    cin >> n_file;
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

    // ******************************************
    // inizializzo
    Real tau = 1;       ///tau iniz.
    Real tau_new = tau;
    Real toll = 1e-12;
    Real r = 0.5;
    Real err = 0;
    int n = 0;  ///conta il numero dei passi
    int f_cnt = 0;  ///conta le chiamate di f

    Real t,T,u[d];   //in u sol. approssimata
    dati(&t,&T,u);


    // scelta del metodo *****************************
    char metodo[]="DP87";// scelta del metodo RK
    const int s=13; // numero degli stadi del metodo scelto;
    const Real p = 8.0; ///ordine del metodo

    // dimensiono array di Butcher *********************************
    Real b[s]={0};
    Real bc[s]={0};
    Real c[s]={0};
    Real A[s][s]={{0}}; //
    // assegno valori array di Butcher
    set_Runge_Kutta_ad(metodo,A[0],b,bc,c);     //bc=b^

    Real diffb[s];
    for (int i=0; i<s; i++) {
        diffb[i] = b[i]-bc[i];
    }

    // inizio il ciclo sul tempo
    stampa_estr(d,t,u,0,&prt);   //stampo valori iniziali
    stampa_estr(d,t,u,0,0);


    while (t<T) {

        Real K[s][d] = {{0}};
        effe(K[0],t,u);
        f_cnt++;

        for (int i=1; i<s; i++) {
            Real V[d] = {0};
            matmat(A[i],K[0],V,1,i,d);

            for (int k=0; k<d; k++) {
                V[k] = u[k] + tau*V[k];
            }
            effe(K[i],t+c[i]*tau,V);        ///Ki finali
            f_cnt++;
        }

        /// stima err
        Real kb[d];
        matmat(diffb,K[0],kb,1,s,d);
        err = tau*(norm2(kb,d));

        tau_new = (pow(r*toll/err, 1.0/p))*tau;

        if (err <= toll) {
            Real aux[d] = {0};
            matmat(b,K[0],aux,1,s,d);   ///F(t,u,tau)

            for (int k=0; k<d; k++) {
                u[k] = u[k] + tau*aux[k];
            }

            t+=tau;
            tau = min(tau_new, T-t);
            n++;

            stampa_estr(d,t,u,err,&prt);
            stampa_estr(d,t,u,err,0);
        }
        else {
            tau = tau_new;
        }
    }

    cout << "cnt passi = " << n << endl;
    cout << "cnt f = " << f_cnt << endl;

    return 0;
}
