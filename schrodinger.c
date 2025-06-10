#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <complex.h>

#define PI 3.14159265358979323846

int main(){
    int N = 1000; // Number of points
    int nciclos = 5; //Numero de ciclos
    double lambda = 1.0; // Wavelength

    double h = 1.0/N;
    double k0_tilde = 2.0 * PI * nciclos / N;
    double s_tilde = 1.0 / (4.0 * k0_tilde * k0_tilde);

    double complex *phi = malloc((N+1) * sizeof(double complex));
    double complex *phi_aux = malloc((N+1) * sizeof(double complex));
    double *V = malloc((N+1) * sizeof(double));

    for(int i=0; i<=N; ++i){
        if(i>=2*N/5 && i <= 3* N /5){
            V[i] = lambda *k0_tilde * k0_tilde;
        }else{
            V[i]=0.0;
        }
    }

    double x0 = N/4.0;
    double sigma = N / 16.0;

    for(int i=0; i<=N; ++i){
        double fase = k0_tilde * i;
        double gauss = exp(-8.0 * pow(4.0 * i - N, 2)/(N *N));
        phi[i] = cexp(I * fase) * gauss;
    }

    phi[0] = 0.0;
    phi[N] = 0.0;

    free(phi);
    free(phi_aux);
    free(V);

    return 0;
}