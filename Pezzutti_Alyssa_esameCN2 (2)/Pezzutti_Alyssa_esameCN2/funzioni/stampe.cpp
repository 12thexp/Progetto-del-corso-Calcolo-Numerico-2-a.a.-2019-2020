#include <iostream>
#include <fstream>
#include "/home/squiddy/Downloads/Calcolo_2_LAB/funzioni/stampe.h"
using namespace std;
typedef double Real;
void stampamat(Real *a,int n,int m,ofstream *prt)
{
    cout.precision(7);
    if (prt==0){
       for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
           cout << *(a+j+i*m) << "\t";
        }
        cout << "\n";
    }
    cout<<" "<<endl;
    }
    else {
        for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
           *prt << *(a+j+i*m) <<"\t";
        }
        *prt << "\n";
    }
    *prt<< " "<<endl;
    }
}

void stampa(int d,Real t,Real *u,ofstream *prt)
{
    if (prt==0){
        //cout << tempo | approx componente 1 | approx componente 2 << endl;
        cout << t << "\t\t";
        for (int i=0;i<d;i++){
            cout << u[i]<<"\t";
        }
                cout <<endl;
    }
    else {
        *prt << t<< "\t";
        for (int i=0;i<d;i++){
            *prt << u[i]<<"\t";
        }
                *prt <<endl;
    }
}


void stampa_estr(int d,Real t,Real *u, Real fc, ofstream *prt)   //fc = fattore di correzione
{
    if (prt==0){
        //cout << tempo | approx componente 1 | approx componente 2 << endl;
        cout << t << "\t\t";
        for(int i=0; i<d; i++){
            cout << u[i] <<"\t";
        }
                cout << fc << endl;
    }
    else {
        *prt << t<< "\t";
        for(int i=0; i<d; i++){
            *prt << u[i] << "\t";
        }
                *prt << fc << endl;
    }
}



