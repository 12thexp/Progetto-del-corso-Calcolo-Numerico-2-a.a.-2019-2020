using namespace std;
typedef  double Real;
void set_Runge_Kutta_esp(char *metodo,Real*A,Real b[],Real c[]);
void set_Runge_Kutta_ad(char *metodo,Real*A,Real b[],Real bc[],Real c[]);
void set_Runge_Kutta_imp(char *metodo,Real*A,Real b[],Real c[]);
void set_Runge_Kutta_diag_imp(char *metodo,Real*A,Real b[],Real c[], Real ind[]);
void set_multi_step(char *metodo,Real ro[],Real sigma[]);
