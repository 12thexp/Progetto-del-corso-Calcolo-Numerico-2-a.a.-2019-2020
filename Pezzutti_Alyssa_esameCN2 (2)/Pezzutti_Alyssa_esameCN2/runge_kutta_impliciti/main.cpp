#include <iostream>
#include <math.h>
#include <fstream>
#include"../../funzioni/mat_vet.h"
#include "../../funzioni/stampe.h"
#include "../../funzioni/problemi.h"
#include "../../funzioni/Runge_Kutta_matrici.h"
#include "../../funzioni/matrici.h"
using namespace std;
typedef double Real;

int main () {

    // predispongo file stampa
    char n_file[21]= {0};
    cout << "dammi nome file di stampa(max 20 caratteri)";
    cin >> n_file;
    ofstream prt(n_file);
    prt.precision(20);
    // definisco puntatori al tipo di funzioni che definiscono il problema
    void (*effe)(Real *,Real,Real *);
    void (*dati)(Real *,Real *, Real *);
    void (*jac)(Real *,Real,Real *);
    //
    // dati del problema **********************************   //le uniche righe che devo cambiare
    const int d=4;
    effe=f_pendolo;
    dati=dati_pendolo;
    jac=jac_pendolo;
    //
    // ******************************************

     // scelta del metodo *****************************
    char metodo[]="Gauss3";// scelta del metodo RK
    const int s=3; // numero degli stadi del metodo scelto;

    // inizializzo
    Real t,T,u[d];   //in u sol. approssimata
    Real delta1, delta2;
    Real delta[d*s];
    Real Jf[d][d]; Real JF[d*s][d*s];
    Real K[s][d] = {{0}};
    Real K2[s][d] = {{0}};
    int f_cnt = 0;

    Real toll = 1e-9;   ///
    int kmax = 200;

    dati(&t,&T,u);

    // dimensiono array di Butcher *********************************
    Real b[s]={0};
    Real c[s]={0};
    Real A[s][s]={{0}}; //
    // assegno valori array di Butcher
    set_Runge_Kutta_imp(metodo,A[0],b,c);

    // N e tau
    unsigned long N;
    cout << "inserisci numero di passi N  = ";
    cin >> N;
    Real tau=(T-t)/Real(N);

    // inizio il ciclo sul tempo
    stampa(d,t,u,&prt);   //stampo valori iniziali
    stampa(d,t,u);

    int i=0, j=0;

    for (i=0; i<s; i++) {
        effe(K[i], t+c[i]*tau, u);
        f_cnt++;
    }


    for (unsigned long n=0; n<=N-1; n++) {
        delta1 = 1e200;
        delta2 = 0;

        for (i=0; i<s; i++) {         ///calcolo la jacobiana
            Real V[d] = {0};
            matmat(A[i],K[0],V,1,s,d);

            for (int j=0; j<d; j++)
                V[j] = u[j]+tau*V[j];

            jac(Jf[0], t+c[i]*tau, V);
            f_cnt+=d;

            for (j=0; j<s; j++) {
                Real aij = -A[i][j]*tau;

                for (int l=0; l<d; l++){
                    for (int m=0; m<d; m++){
                        JF[i*d+l][j*d+m] = aij*Jf[l][m];
                    }
                }
            }
        }
        for (int k=0; k<d*s; k++)
            JF[k][k]+=1;

        int P[d*s] = {0};
        lu(JF[0],P,d*s);       ///fattorizzo JF

        int k=0;
        //int newt = 0;
        do {
            for (i=0; i<s; i++) {
                Real W[d] = {0};
                matmat(A[i],K[0],W,1,s,d);

                for (j=0; j<d; j++) {
                    W[j] = u[j]+tau*W[j];
                }

                effe(K2[i], t+c[i]*tau, W);
                f_cnt++;
            }
            int iKK=0;
            Real tn[d*s];
            for (int l=0; l<s; l++){       ///costruisco tn = K2-K
                for (int m=0; m<d; m++){
                    tn[iKK] = K2[l][m] - K[l][m];
                    iKK++;
                }
            }
            risist(JF[0],P,delta,tn,d*s);            ///dim: JF=ds*ds, delta=ds*1,tn=ds*1

            iKK=0;              ///K+= delta
            for (int l=0; l<s; l++) {
                for (int m=0; m<d; m++) {
                    K[l][m] += delta[iKK];
                    iKK++;
                }
            }
            k++;
            delta2 = norm2(delta,d*s);
            if (delta2 > delta1)
                break;
            else
                delta1 = delta2;
        //newt++;
        }
        while ((delta1>toll) && (k<kmax));
//        cout << t << '\t';
//        cout << newt << endl;

        Real aux[d] = {0};
        matmat(b,K[0],aux,1,s,d);

        for (int kk=0; kk<d; kk++) {
            u[kk] = u[kk] + tau*aux[kk];
        }
        t+=tau;
        stampa(d,t,u,&prt);
        stampa(d,t,u,0);
    }
    cout << "f_cnt:" << f_cnt << endl;
    return 0;
}





