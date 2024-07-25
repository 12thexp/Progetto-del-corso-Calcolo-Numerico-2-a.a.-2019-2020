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
    prt.precision(30);
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
    Real t,T,u[d];   //in u sol. approssimata
    dati(&t,&T,u);
    int f_cnt = 0;  ///conta le chiamate di f

    // scelta del metodo *****************************
    char metodo[]="RK4";// scelta del metodo RK
    const int s=4; // enumero degli stadi del metodo scelto;

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
    stampa(d,t,u,&prt);   //stampo valori iniziali
    stampa(d,t,u);

    for (unsigned long n=1; n<=N; n++) {
        Real K[s][d] = {{0}};
        effe(K[0],t,u);
        f_cnt++;


        for (int i=1; i<s; i++) {
            Real V[d] = {0};
            matmat(A[i],K[0],V,1,i,d);

            for (int kk=0; kk<d; kk++) {
                V[kk] = u[kk] + tau*V[kk];
            }
            effe(K[i],t+c[i]*tau,V);
            f_cnt++;
        }

        Real aux[d] = {0};
        matmat(b,K[0],aux,1,s,d);

        for (int kk=0; kk<d; kk++) {
            u[kk] = u[kk] + tau*aux[kk];
        }

        t+=tau;
//
//         if (n == 1) {
//            stampa(d,t,u,&prt);     //stampa primo passo
//            stampa(d,t,u);
//        }

        stampa(d,t,u,&prt);
        stampa(d,t,u);
    }
//    stampa(d,t,u,&prt);     //stampa finale
//    stampa(d,t,u);
    cout << "cnt f = " << f_cnt << endl;

    return 0;
}
