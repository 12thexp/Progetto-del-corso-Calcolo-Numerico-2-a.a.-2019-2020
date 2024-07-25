using namespace std;
typedef double Real;
#include <iostream>
#include <string.h>
#include <math.h>
#include "Runge_Kutta_matrici.h"
void set_Runge_Kutta_esp(char *metodo,Real*A,Real b[],Real c[]) {
    int s;
    int m=0;
    if(strcmp(metodo,"EE")==0){m=1;s=1;}  ///m = numero metodo per switch case
    if(strcmp(metodo,"Heun")==0){m=2;s=2;}
    if(strcmp(metodo,"RK4")==0){m=3;s=4;}
    if(strcmp(metodo,"es6_3")==0){m=4;s=3;}
    if(strcmp(metodo,"Fehlberg")==0){m=5;s=6;}

// assegno valori matrice e vettori
        switch(m) {
            case 1:{
    //          Eulero Esplicito
                b[0]=1.0;
            }
            break;
            case 2:{
    //          Heun
                b[0]=b[1]=0.5;
                *(A+2)=1.0;
                }
            break;
            case 3:{
    //          RK4
                b[0]=b[3]=1.0/6.0;
                b[1]=b[2]=1.0/3.0;
                *(A+4)=*(A+9)=0.5;
                *(A+14)=1.0;
                }
            break;
            case 4:{    //gli elementi di A si contano da 0!!
    //          es6_3
                b[0]=2.0/9.0;
                b[1]=1.0/3.0;
                b[2]=4.0/9.0;
                *(A+3)=0.5;
                *(A+7)=3.0/4.0;
                }
            break;
            case 5:{
        //      Fehlberg
                b[0]=16.0/135.0;
                b[1]=0.0;
                b[2]=6656.0/12825.0;
                b[3]=28561.0/56430.0;
                b[4]=-9.0/50.0;
                b[5]=2.0/55.0;
                //riga 2
                *(A+6)=0.25;
                //riga 3
                *(A+12)=3.0/32.0;
                *(A+13)=9.0/32.0;
                //riga 4
                *(A+18)=1932.0/2197.0;
                *(A+19)=-7200.0/2197.0;
                *(A+20)=7296.0/2197.0;
                //riga 5
                *(A+24)=439.0/216.0;
                *(A+25)=-8.0;
                *(A+26)=3680.0/513.0;
                *(A+27)=-845.0/4104.0;
                //riga 6
                *(A+30)=-8.0/27.0;
                *(A+31)=2.0;
                *(A+32)=-3544.0/2565.0;
                *(A+33)=1859.0/4104.0;
                *(A+34)=-11.0/40.0;
                }
            break;
            }
// calcolo C
        for(int k=0;k<s;k++){
            Real sum=*(A+k*s);
            for (int j=1;j<s;j++){
                sum+=*(A+s*k+j);
            }
            c[k]=sum;
         }
}


///
void set_Runge_Kutta_imp(char *metodo,Real*A,Real b[],Real c[]){
    int ns;
    int m=0;
    if(strcmp(metodo,"Gauss2")==0){m=2;ns=2;}
    if(strcmp(metodo,"Gauss3")==0){m=3;ns=3;}
    if(strcmp(metodo,"EI")==0){m=4;ns=1;}
    if(strcmp(metodo,"Gauss1")==0){m=1;ns=1;}
    if(strcmp(metodo,"Trapezi")==0){m=5;ns=2;}

// inizializzo a zero matrice A e il vettore Ind a 1

// assegno valori matrice e vettori
    switch(m) {
        case 1:{
    // Gauss 1
            b[0]=1;
            *A=0.5;
            }
        break;
        case 2:{
    // Gauss2
            b[0]=0.5;
            b[1]=0.5;
            *A=1.0/4.0;
            *(A+1)=(3.0-2.0*sqrt(3))/12.0;
            *(A+2)=(3.0+2.0*sqrt(3))/12.0;
            *(A+3)=1.0/4.0;
            }
        break;
        case 3:{
    // Gauss3
            b[0]=5.0/18.0;
            b[1]=4.0/9.0;
            b[2]=5.0/18.0;
            *A=5.0/36.0;
            *(A+1)=2.0/9.0-sqrt(15)/15.0;
            *(A+2)=5.0/36.0-sqrt(15)/30.0;
            *(A+3)=5.0/36.0+sqrt(15)/24.0;
            *(A+4)=2.0/9.0;
            *(A+5)=5.0/36.0-sqrt(15)/24.0;
            *(A+6)=5.0/36.0+sqrt(15)/30.0;
            *(A+7)=2.0/9.0+sqrt(15)/15.0;
            *(A+8)=5.0/36.0;
            }
        break;
        case 4:{
    // Eulero Implicito
            b[0]=1;
            *A=1.0;
            }
        break;
        case 5:{
    // Gauss2
            b[0]=0.5;
            b[1]=0.5;
            *A=0.0;
            *(A+1)=0.0;
            *(A+2)=1.0/2.0;
            *(A+3)=1.0/2.0;
            }
        break;
    }
       // calcolo C
          for(int k=0;k<ns;k++){
            Real sum=*(A+k*ns);
            for (int j=1;j<ns;j++){
            sum+=*(A+ns*k+j);
            }
        c[k]=sum;
         }
}



///
void set_Runge_Kutta_diag_imp(char *metodo,Real*A,Real b[],Real c[], Real ind[]){
    int ns;
    int m=0;
    if(strcmp(metodo,"EI")==0){m=1;ns=1;}
    if(strcmp(metodo,"Trap_imp")==0){m=2;ns=2;}

// inizializzo a zero matrice A e il vettore Ind a 1

// assegno valori matrice e vettori
    switch(m) {
        case 1:{
    // Eulero Implicito
            b[0]=1;
            *A=1.0;
            }
        break;
        case 2:{
    // Trapezi impliciti
            b[0]=0.5;
            b[1]=0.5;
            *A=0.0;
            *(A+1)=0.0;
            *(A+2)=1.0/2.0;
            *(A+3)=1.0/2.0;
            }
        break;

        }
       // calcolo C
    for(int k=0;k<ns;k++) {
        Real sum=*(A+k*ns);
        for (int j=1;j<ns;j++){
            sum+=*(A+ns*k+j);
        }
        c[k]=sum;
    }


    for (int k=0; k<ns*ns; k++) {
        if (*(A+k) != 0.0)
            ind[k] = 1;
        }
    }




///
void set_Runge_Kutta_ad(char *metodo,Real*A,Real b[],Real bc[],Real c[]) {       //bc = b^, coeff del metodo p-1
    int ns;
    int m=0;
    if(strcmp(metodo,"DP87")==0){m=1;ns=13;}
    if(strcmp(metodo,"RKF54")==0){m=2;ns=6;}
    if(strcmp(metodo,"EEHEUN")==0){m=3;ns=2;}
// assegno valori matrice e vettori
        switch(m) {
//    // **********************************************************
    case 1:{
//        DP87
        for(int j=0;j<ns;j++){
        bc[j]=b[j]=0.0;
        }
            b[0]=14005451.0/335480064.0;
            b[5]=-59238493.0/1068277825.0;
            b[6]=181606767.0/758867731.0;
            b[7]=561292985.0/797845732.0;
            b[8]=-1041891430.0/1371343529.0;
            b[9]=760417239.0/1151165299.0;
            b[10]=118820643.0/751138087.0;
            b[11]=-528747749.0/2220607170.0;
            b[12]=0.25;
            bc[0]=13451932.0/455176623.0;
            bc[5]=-808719846.0/976000145.0;
            bc[6]=1757004468.0/5645159321.0;
            bc[7]=656045339.0/265891186.0;
            bc[8]=-3867574721.0/1518517206.0;
            bc[9]=465885868.0/322736535.0;
            bc[10]=53011238.0/667516719.0;
            bc[11]=2.0/45.0;
            c[0]=0.0;
            c[1]=1.0/18.0;
            c[2]=1.0/12.0;
            c[3]=0.125;
            c[4]=5.0/16.0;
            c[5]=3.0/8.0;
            c[6]=59.0/400.0;
            c[7]=93.0/200.0;
            c[8]=5490023248.0/9719169821.0;
            c[9]=13.0/20.0;
            c[10]=1201146811.0/1299019798.0;
            c[11]=1.0;
            c[12]=1.0;
        //riga 2
        *(A+13)=1.0/18.0;
        //riga 3
        *(A+26)=1.0/48.0;
        *(A+27)=1.0/16.0;
        //riga 4
        *(A+39)=1.0/32.0;
        *(A+41)=3.0/32.0;
        //riga 5
        *(A+52)=5.0/16.0;
        *(A+54)=-75.0/64.0;
        *(A+55)=75.0/64.0;
        //riga 6
        *(A+65)=3.0/80.0;
        *(A+68)=3.0/16.0;
        *(A+69)=3.0/20.0;
        //riga 7
        *(A+78)=29443841.0/614563906.0;
        *(A+81)=77736538.0/692538347.0;
        *(A+82)=-28693883.0/1125000000.0;
        *(A+83)=23124283.0/1800000000.0;
        //riga 8
        *(A+91)=16016141.0/946692911.0;
        *(A+94)=61564180.0/158732637.0;
        *(A+95)=22789713.0/633445777.0;
        *(A+96)=545815736.0/2771057229.0;
        *(A+97)=-180193667.0/1043307555.0;
        //riga 9
        *(A+104)=39632708.0/573591083.0;
        *(A+107)=-433636366.0/683701615.0;
        *(A+108)=-421739975.0/2616292301.0;
        *(A+109)=100302831.0/723423059.0;
        *(A+110)=790204164.0/839813087.0;
        *(A+111)=800635310.0/3783071287.0;
        //riga 10
        *(A+117)=246121993.0/1340847787.0;
        *(A+120)=-37695042795.0/15268766246.0;
        *(A+121)=-309121744.0/1061227803.0;
        *(A+122)=-12992083.0/490766935.0;
        *(A+123)=6005943493.0/2108947869.0;
        *(A+124)=393006217.0/1396673457.0;
        *(A+125)=123872331.0/1001029789.0;
        //riga 11
        *(A+130)=-1028468189.0/846180014.0;
        *(A+133)=8478235783.0/508512852.0;
        *(A+134)=1311729495.0/1432422823.0;
        *(A+135)=-10304129995.0/1701304382.0;
        *(A+136)=-48777925059.0/3047939560.0;
        *(A+137)=15336726248.0/1032824649.0;
        *(A+138)=-45442868181.0/3398467696.0;
        *(A+139)=3065993473.0/597172653.0;
        //riga 12
        *(A+143)=185892177.0/718116043.0;
        *(A+146)=-3185094517.0/667107341.0;
        *(A+147)=-477755414.0/1098053517.0;
        *(A+148)=-703635378.0/230739211.0;
        *(A+149)=5731566787.0/1027545527.0;
        *(A+150)=5232866602.0/850066563.0;
        *(A+151)=-4093664535.0/808688257.0;
        *(A+152)=3962137247.0/1805957418.0;
        *(A+153)=65686358.0/487910083.0;
        //riga 13
       *(A+156)=403863854.0/491063109.0;
        *(A+159)=-5068492393.0/434740067.0;
        *(A+160)=-411421997.0/543043805.0;
        *(A+161)=652783627.0/914296604.0;
        *(A+162)=11173962825.0/925320556.0;
        *(A+163)=-13158990841.0/6184727034.0;
        *(A+164)=3936647629.0/1978049680.0;
        *(A+165)=-160528059.0/685178525.0;
        *(A+166)=248638103.0/1413531060.0;
        }
    break;
    // **********************************************************
    case 2:{
//      RKF54
        b[0]=16.0/135.0;
        b[1]=0.0;
        b[2]=6656.0/12825.0;
        b[3]=28561.0/56430.0;
        b[4]=-9.0/50.0;
        b[5]=2.0/55.0;
        bc[0]=25.0/216.0;
        bc[1]=0.0;
        bc[2]=1408.0/2565.0;
        bc[3]=2197.0/4104.0;
        bc[4]=-1.0/5.0;
        bc[5]=0.0;
        //riga 2
        *(A+6)=0.25;
        //riga 3
        *(A+12)=3.0/32.0;
        *(A+13)=9.0/32.0;
        //riga 4
        *(A+18)=1932.0/2197.0;
        *(A+19)=-7200.0/2197.0;
        *(A+20)=7296.0/2197.0;
        //riga 5
        *(A+24)=439.0/216.0;
        *(A+25)=-8.0;
        *(A+26)=3680.0/513.0;
        *(A+27)=-845.0/4104.0;
        //riga 6
        *(A+30)=-8.0/27.0;
        *(A+31)=2.0;
        *(A+32)=-3544.0/2565.0;
        *(A+33)=1859.0/4104.0;
        *(A+34)=-11.0/40.0;
        }
    break;

    case 3:{
//       EEHEUN
        b[0]=0.5;
        b[1]=0.5;
        bc[0]=1.0;
        bc[1]=0.0;
        //riga 2
        *(A+2)=1.0;
        }
    break;
    }
    // **********************************************************
        // calcolo C
          for(int k=0;k<ns;k++){
            Real sum=*(A+k*ns);
            for (int j=1;j<ns;j++){
            sum+=*(A+ns*k+j);
            }
        c[k]=sum;
  //      cout << c[k]<<endl;
         }
}



///multistep
using namespace std;
typedef double Real;
#include <iostream>
#include <string.h>

//
void set_multi_step(char *metodo,Real ro[],Real sigma[]){
    int m=0;
    // AB sta per Adams_Bashforth, AM per Adams_Moulton
    // il numero che segue corrisponde all'ordine del metodo
    if(strcmp(metodo,"AB2")==0)m=1;//passi=2
    if(strcmp(metodo,"AB3")==0)m=2;//passi=3;}
    if(strcmp(metodo,"AB4")==0)m=3;//passi=4;}
    if(strcmp(metodo,"AB1")==0)m=8;//passi=1;}
    if(strcmp(metodo,"AM1")==0)m=7;//passi=1;}
    if(strcmp(metodo,"AM2")==0)m=4;//passi=1;}
    if(strcmp(metodo,"AM3")==0)m=5;//passi=2;}
    if(strcmp(metodo,"AM4")==0)m=6;//passi=3;}
    if(strcmp(metodo,"BDF2")==0)m=9;//passi=2;}
    if(strcmp(metodo,"BDF3")==0)m=10;//passi=3;}
    if(strcmp(metodo,"BDF4")==0)m=11;//passi=4;}
    if(strcmp(metodo,"BDF5")==0)m=13;//passi=5;}
    if(strcmp(metodo,"BDF6")==0)m=14;//passi=6;}
    if(strcmp(metodo,"LF")==0)m=12;//passi=2;}
    //

    switch(m){
            case 0:{
                cout <<metodo<< "  metodo non previsto abortisci processo"<<endl;
                }
                break;
            case 1:{
                ro[0]=1.0;
                ro[1]=-1.0;
                ro[2]=0.0;
                sigma[0]=0.0;
                sigma[1]=3.0/2.0;
                sigma[2]=-0.5;
                }
                break;
            case 2:{
                ro[0]=12.0;
                ro[1]=-12.0;
                ro[2]=0.0;
                ro[3]=0.0;
                sigma[0]=0.0;
                sigma[1]=23.0;
                sigma[2]=-16.0;
                sigma[3]=5.0;
                }
                break;
            case 3:{
                ro[0]=24.0;
                ro[1]=-24.0;
                ro[2]=0.0;
                ro[3]=0.0;
                ro[4]=0.0;
                sigma[0]=0.0;
                sigma[1]=55.0;
                sigma[2]=-59.0;
                sigma[3]=37.0;
                sigma[4]=-9.0;
                }
                break;
            case 4:{
                ro[0]=1.0;
                ro[1]=-1.0;
                sigma[0]=0.5;
                sigma[1]=0.5;
                }
                break;
            case 5:{
                ro[0]=12.0;
                ro[1]=-12.0;
                ro[2]=0.0;
                sigma[0]=5.0;
                sigma[1]=8.0;
                sigma[2]=-1.0;
                }
                break;
            case 6:{
                ro[0]=24.0;
                ro[1]=-24.0;
                ro[2]=0.0;
                ro[3]=0.0;
                sigma[0]=9.0;
                sigma[1]=19.0;
                sigma[2]=-5.0;
                sigma[3]=1.0;
                }
                break;
            case 7:{
                ro[0]=1.0;
                ro[1]=-1.0;
                sigma[0]=1.0;
                sigma[1]=0.0;
                }
                break;
            case 8:{
                ro[0]=1.0;
                ro[1]=-1.0;
                sigma[0]=0.0;
                sigma[1]=1.0;
                }
                break;
            case 9:{
                ro[0]=1.5;
                ro[1]=-2.0;
                ro[2]=0.5;
                sigma[0]=1.0;
                sigma[1]=0.0;
                sigma[2]=0.0;
                }
                break;
            case 10:{
                ro[0]=11.0/6.0;
                ro[1]=-3.0;
                ro[2]=1.5;
                ro[3]=-1.0/3.0;
                sigma[0]=1.0;
                sigma[1]=0.0;
                sigma[2]=0.0;
                sigma[3]=0.0;
                }
                break;
            case 11:{
                ro[0]=25.0/12.0;
                ro[1]=-4.0;
                ro[2]=3.0;
                ro[3]=-4.0/3.0;
                ro[4]=0.25;
                sigma[0]=1.0;
                sigma[1]=0.0;
                sigma[2]=0.0;
                sigma[3]=0.0;
                sigma[4]=0.0;
                }
                break;
            case 12:{
                ro[0]=1.0;
                ro[1]=0.0;
                ro[2]=-1.0;
                sigma[0]=0.0;
                sigma[1]=2.0;
                sigma[2]=0.0;
                }
                break;
            case 13:{
                ro[0]=137.0/60.0;
                ro[1]=-5.0;
                ro[2]=5.0;
                ro[3]=-10.0/3.0;
                ro[4]=1.25;
                ro[5]=-0.2;
                sigma[0]=1.0;
                sigma[1]=0.0;
                sigma[2]=0.0;
                sigma[3]=0.0;
                sigma[4]=0.0;
                sigma[5]=0.0;
                }
                break;
            case 14:   {
				ro[0]=49.0/20.0;
                ro[1]=-6.0;
                ro[2]=15.0/2.0;
                ro[3]=-20.0/3.0;
                ro[4]=15.0/4.0;
                ro[5]=-6.0/5.0;
                ro[6]=1.0/6.0;
                sigma[0]=1.0;
                sigma[1]=0.0;
                sigma[2]=0.0;
                sigma[3]=0.0;
                sigma[4]=0.0;
                sigma[5]=0.0;
                sigma[6]=0.0;
                }
                break;
    }
}

