#include "matrici.h"
#include <math.h>
using namespace std;
typedef  double Real;

void vander(Real *A,int n){
    // matrice di van der Monde costruita sul vettore[1,2,...n]
    int k=0;
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
           A[k]=pow(i+1.0,j);
           k++;
        }
    }
        }

void hilbert(Real *A, int n){
    // matrice di hilbert
    int k=0;
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
           A[k]=1.0/(i+j+1.0);
           k++;
        }
    }
}

void es8_1(Real *A, int n) {
    for (int i=0; i<n; i++)
        A[i]=-20;
}
