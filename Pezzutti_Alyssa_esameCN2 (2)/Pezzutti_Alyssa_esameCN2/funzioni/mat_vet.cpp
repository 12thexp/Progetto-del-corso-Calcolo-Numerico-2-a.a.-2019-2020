#include "mat_vet.h"
#include <iostream>
#include <math.h>
using namespace std;
typedef  double Real;

void gauss(Real *A, Real x[], Real b[], int n)

// risolve sistema lineare con eliminazione gaussiana con pivot parziale
// n dimensione corrente della matrice e di b
{
    Real aux;
//
    for (int k=0;k<n-1;k++){
//  cerco il pivot sulla colonna
        Real pivot=fabs(*(A+k*n+k));
        int ipiv=k;
        for(int kk=k+1;kk<n;kk++){
            aux=fabs(*(A+kk*n+k));
            if(aux>pivot){
                ipiv=kk;
                pivot=aux;
            }
        }
//  scambio le righe
        if(ipiv != k){
            for(int kk=k;kk<n;kk++){
                aux=*(A+k*n+kk);
                *(A+k*n+kk)=*(A+ipiv*n+kk);
                *(A+ipiv*n+kk)=aux;
            }
            aux=b[k];
            b[k]=b[ipiv];
            b[ipiv]=aux;
        }
            for (int i=k+1; i<n; i++){
                Real m=*(A+i*n+k)/(*(A+k*n+k));
                *(A+i*n+k)=0.0;
                for(int j=k+1;j<n;j++){
                *(A+i*n+j)=*(A+i*n+j)-m*(*(A+k*n+j));
                }
                b[i]=b[i]-m*b[k];
            }
    }
    x[n-1]=b[n-1]/(*(A+(n-1)*n+n-1));

    for(int i=n-2;i>=0;i--){
        Real sum=0.0;
        for(int j=i+1;j<n;j++){
            sum += *(A+i*n+j)*x[j];
        }
        x[i]=(b[i]-sum)/(*(A+i*n+i));
    }
}



void lu(Real *A, int P[], int n)
// risolve sistema lineare con eliminazione gaussiana con pivot parziale
// n dimensione corrente della matrice e di b
{
    Real aux;
    int ia;
    for (int k=0;k<n;k++){
        P[k]=k;
    }
    for (int k=0;k<n-1;k++){
//  cerco il pivot sulla colonna
        Real pivot=fabs(*(A+k*n+k));
        int ipiv=k;
        for(int kk=k+1;kk<n;kk++){
            aux=fabs(*(A+kk*n+k));
            if(aux>pivot){
                ipiv=kk;
                pivot=aux;
            }
        }
//  scambio le righe
    if(ipiv != k){
        for(int kk=0;kk<n;kk++){
            aux=*(A+k*n+kk);
            *(A+k*n+kk)=*(A+ipiv*n+kk);
            *(A+ipiv*n+kk)=aux;
        }
        ia=P[k];
        P[k]=P[ipiv];
        P[ipiv]=ia;
    }
        for (int i=k+1; i<n; i++){
            Real m=*(A+i*n+k)/(*(A+k*n+k));
            *(A+i*n+k)=m;
            for (int j=k+1;j<n;j++){
            *(A+i*n+j)-=m*(*(A+k*n+j));
            }
        }
    }
}

void risist(Real *A, int P[],Real x[], Real b[], int n){
    Real sum=0.0;
    for (int k=0;k<n;k++){
       x[k]=b[P[k]];
    }
    for (int k=1;k<n;k++){
        sum=0.0;
        for(int i=0;i<k;i++){
        sum+= *(A+k*n+i)*x[i];
        }
        x[k]= x[k]-sum;
    }
    x[n-1]=x[n-1]/(*(A+(n-1)*n+n-1));
    for(int i=n-2;i>=0;i--){
        sum=0.0;
        for(int j=i+1;j<n;j++){
            sum += *(A+i*n+j)*x[j];
        }
        x[i]=(x[i]-sum)/(*(A+i*n+i));
    }
}

void matmat(Real *A, Real *B, Real *Ris, int ra, int ca, int cb)
// calcola prodotto matrice raxca per matrice rbxcb con ca=rb
// A*B = Ris
{
    for (int i=0; i<ra; i++){
        for (int j=0; j<cb; j++){
                Real sum=0.0;
                for (int k=0; k<ca; k++){
                        sum += (*(A+i*ca+k))*(*(B+j+k*cb));
                }
                *(Ris+i*cb+j)=sum;
        }
    }
}


///norme

Real norm2(Real v[], int n) {
    Real N = 0;
    for (int i=0; i<n; i++) {
        N += pow(v[i], 2.0);
    }
    N = sqrt(N);
    return N;
}

///copia vettori
void vcopy(Real from[], Real to[], int l) {
    for (int i=0; i<l; i++)
        to[i] = from[i];
}


