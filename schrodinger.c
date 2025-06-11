#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <complex.h>

#define PI 3.14159265358979323846
void funcionOndaInicial(int N, double h, double x0, double sigma, double k0, double complex *phi);
void tridiagonal(int N, double s, double *V, double complex *a, double complex *b, double complex *c);
void resolverTridiagonal(int N, double complex *a, double complex *b, double complex *c, double complex *d, double complex *x);
void evolucionarOnda(int N, int pasos, double s, double *V, double complex *phi);

int main(){
    int N = 1000; // Number of points
    int nciclos = 5; //Numero de ciclos
    double lambda = 1.0; // Wavelength

    double h = 1.0/N;
    double k0 = 2.0*PI*nciclos/(N*h);

    double complex *phi = malloc((N+1) * sizeof(double complex));
    double complex *phi_aux = malloc((N+1) * sizeof(double complex));
    double *V = malloc((N+1) * sizeof(double));

    for(int i=0; i<=N; ++i){
        if(i>=2*N/5 && i <= 3* N /5){
            V[i] = lambda *k0 * k0;
        }else{
            V[i]=0.0;
        }
    }

    double x0=N/4.0*h;
    double sigma=N/16.0*h;
    funcionOndaInicial(N, h, x0, sigma, k0, phi);

    free(phi);
    free(phi_aux);
    free(V);

    return 0;
}

void funcionOndaInicial(int N, double h, double x0, double sigma, double k0, double complex *phi){
    double norma=0.0;
    
    for(int i=0; i<=N; ++i){
        double x=i*h;
        double gauss = exp(-0.5*pow((x - x0) / sigma, 2));
        double complex fase=cexp(I*k0*x);
        phi[i]=gauss*fase;
        norma+= creal(phi[i]*conj(phi[i]));
    }

    norma = sqrt(norma);
    for(int i=0; i<=N; ++i){
        phi[i] /= norma;
    }

    phi[0] = 0.0; // Boundary condition
    phi[N] = 0.0; // Boundary condition
}

void tridiagonal(int N, double s, double *V, double complex *a, double complex *b, double complex *c){
    for(int j=1; j<N; ++j){
        a[j-1]=-I*s;
        b[j-1]=1.0 + 2.0*I*s + I*s*V[j];
        c[j-1]=-I*s;
    }
}

void resolverTridiagonal(int N, double complex *a, double complex *b, double complex *c, double complex *d, double complex *x){
    double complex *c_prime = malloc((N-1) * sizeof(double complex));
    double complex *d_prime = malloc(N * sizeof(double complex));

    c_prime[0]= c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    for(int i=1; i<N-1; ++i){
        double complex m = b[i]-a[i-1]*c_prime[i-1];
        c_prime[i] = c[i] / m;
        d_prime[i] = (d[i] - a[i-1]*d_prime[i-1]) / m;
    }

    d_prime[N-1] = (d[N-1] - a[N-2]*d_prime[N-2]) / (b[N-1] - a[N-2]*c_prime[N-2]);

    x[N-1] = d_prime[N-1];
    for(int i=N-2; i>=0; --i){
        x[i] = d_prime[i] - c_prime[i]*x[i+1];
    }

    free(c_prime);
    free(d_prime);
}

void evolucionarOnda(int N, int pasos, double s, double *V, double complex *phi){
    double complex *a = malloc((N-1) * sizeof(double complex));
    double complex *b = malloc((N-1) * sizeof(double complex));
    double complex *c = malloc((N-1) * sizeof(double complex));
    double complex *d = malloc((N-1) * sizeof(double complex));
    double complex *sol = malloc((N-1) * sizeof(double complex));

    tridiagonal(N, s, V, a, b, c);

    for(int t=0; t<pasos; ++t){
        for(int j=1; j<N; ++j){
            d[j-1]=(1.0 - 2.0*I*s - I*s*V[j]) *phi[j] + I*s*(phi[j-1] + phi[j+1]);
        }

        resolverTridiagonal(N-1, a, b, c, d, sol);

        phi[0]=0.0;
        for(int j=1; j<N; ++j){
            phi[j] = sol[j-1];
        }
        phi[N] = 0.0; // Boundary condition
    }

    free(a);
    free(b);
    free(c);
    free(d);
    free(sol);
}

double probabilidadTransmision(int N, double complex *phi){
    int i0 = 3*N/5;
    double suma =0.0;
    for(int i=i0; i<=N; ++i){
        suma += creal(phi[i]*conj(phi[i]));
    }
    return suma;
}

