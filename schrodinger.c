#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define PI 3.14159265358979323846
#define N 500
#define h 0.01
#define T 3000
#define M 100
typedef double complex cplx;

void inicializarPhi(cplx *phi, int nciclos);
void inicializarPotencial(double *V, double lambda, double k_tilde);
void construirTridiagonal(cplx *A1, cplx *A2, cplx *A3, double *V, double s_tilde);
void alfaGamma(cplx *A1, cplx *A2, cplx *A3, cplx *alfa, cplx *gamma);
void pasoTemporal(cplx *phi, cplx *chi, cplx *A3, cplx *b, cplx *gamma, cplx *alfa, cplx *beta, double s_tilde, cplx *phi_next);
double probabilidad(cplx *phi);
double norma(cplx *phi);
int maximo(double *PD, int t);
double valor_esperado_x(cplx *phi);
double valor_esperado_T(cplx *phi);
void guardarArray(const char *nombre, double *array, int n);
void guardarPerfil(const char *nombre, cplx *phi, double *V);

void inicializarPhi(cplx *phi, int nciclos){
    double k_tilde = 2.0*PI*nciclos/N;
    double x0 = N*h/4.0;
    double sigma = N*h/16.0;
    double norm = 0.0;

    for(int j=0; j<=N; ++j){
        double x = j*h;
        phi[j] = cexp(I*k_tilde*x)*exp(-0.5 * pow((x-x0)/sigma, 2));
        norm += pow(cabs(phi[j]), 2)*h;
    }

    norm = sqrt(norm);
    for(int j=0; j<=N; ++j){
        phi[j] /= norm;
    }

    phi[0] = 0.0; // Boundary condition
    phi[N] = 0.0; // Boundary condition
}

void inicializarPotencial(double *V, double lambda, double k_tilde){
    int inicio = 2*N/5;
    int fin = 3*N/5;
    double V0 = lambda *k_tilde*k_tilde;

    for(int j=0; j<=N; ++j){
        if(j>=inicio && j<=fin){
            V[j] = V0;
        } else {
            V[j] = 0.0;
        }
    }
}

void construirTridiagonal(cplx *A1, cplx *A2, cplx *A3, double *V, double s_tilde){
    for(int j=1; j<N; ++j){
        A1[j] = 1.0;
        A2[j] = -2.0 + 2.0*I/s_tilde - V[j];
        A3[j] = 1.0;
    }
}

void alfaGamma(cplx *A1, cplx *A2, cplx *A3, cplx *alfa, cplx *gamma){
    alfa[N-1] = 0.0;
    for(int j=N-1; j>0; --j){
        gamma[j] = A2[j] + A3[j] *alfa[j];
        alfa[j-1] = -A1[j]/gamma[j];
    }
}

void pasoTemporal(cplx *phi, cplx *chi, cplx *A3, cplx *b, cplx *gamma, cplx *alfa, cplx *beta, double s_tilde, cplx *phi_next){
    for(int j=1; j<N; ++j){
        b[j] = 4.0*I/s_tilde*phi[j];
    }

    beta[N-1] = 0.0;
    for(int j=N-1; j>0; --j){
        beta[j-1] = (b[j]-A3[j]*beta[j])/gamma[j];
    }

    chi[0] = 0.0;
    chi[N] = 0.0;
    for(int j=0; j<N; ++j){
        chi[j+1] = alfa[j]*chi[j] + beta[j];
    }

    for(int j=1; j<N; ++j){
        phi_next[j] = chi[j] - phi[j];
    }

    for(int j=1; j<N; ++j){
        phi[j] = phi_next[j];
    }
}

double probabilidad(cplx *phi){
    double suma = 0.0;
    int ini = 3*N/5 + 1; // Justo después de la barrera
    for(int j=ini; j<=N; ++j){
        suma += pow(cabs(phi[j]), 2)*h;
    }
    return suma;
}

double norm(cplx *phi){
    double norma = 0.0;
    for(int j=0; j<=N; ++j){
        norma += pow(cabs(phi[j]), 2)*h;
    }
    return norma;
}

int maximo(double *PD, int t){
    for(int n=1; n<t; ++n){
        if(PD[n] > PD[n-1] && PD[n] > PD[n+1]){
            return n;
        }
    }
    return -1; // No maximum found
}

double valor_esperado_x(cplx *phi) {
    double suma = 0.0;
    for (int j = 0; j <= N; j++) {
        double x = j * h;
        suma += x * pow(cabs(phi[j]), 2) * h;
    }
    return suma;
}

double valor_esperado_T(cplx *phi) {
    double suma = 0.0;
    for (int j = 1; j < N; j++) {
        cplx laplaciano = (phi[j+1] - 2.0*phi[j] + phi[j-1]) / (h*h);
        suma += creal(conj(phi[j]) * laplaciano) * h;
    }
    return -0.5 * suma;
}


void guardarArray(const char *nombre, double *array, int n){
    FILE *f = fopen(nombre, "w");
    if (f == NULL) {
        perror("Error opening file");
        return;
    }

    for(int i=0; i<n; ++i){
        fprintf(f, "%d %.10f\n", i, array[i]);
    }
    fclose(f);  
}

void guardarPerfil(const char *nombre, cplx *phi, double *V){
    FILE *f = fopen(nombre, "w");
    if (f == NULL) {
        perror("Error opening file");
        return;
    }
    for(int j=0; j<=N; ++j){
        double x = j * h;
        double densidad = pow(cabs(phi[j]), 2);
        fprintf(f, "%.8f %.10f %.10f\n", x, densidad, V[j]);
    }
    fclose(f);
}

int main() {
    int nciclos = 20;
    double lambda_min = 0.1;
    double lambda_max = 2.0;
    double lambda_step = 0.05;

    const char *output_dir = "C:\\Users\\molin\\Escritorio\\Universidad\\3Fisica\\Segundocuatri\\Compu\\Schrodinger\\archivos\\";

    for (double lambda = lambda_min; lambda <= lambda_max + 1e-8; lambda += lambda_step) {
        double k_tilde = 2.0 * PI * nciclos / N;
        double s_tilde = 1.0 / (4.0 * k_tilde * k_tilde);

        cplx phi[N+1], chi[N+1], phi_next[N+1];
        cplx A1[N+1], A2[N+1], A3[N+1];
        cplx gamma[N], alfa[N], beta[N], b[N+1];
        double V[N+1];
        double PD[T], normas[T], xesp[T], Tesp[T];

        int mT = 0;
        double PD_total = 0.0;
        srand(time(NULL));

        // ----- Nombres de archivos dinámicos -----
        char nombre_PD[256], nombre_norma[256], nombre_xesp[256], nombre_Tesp[256], nombre_phiV[256], nombre_PDnd[256];
        sprintf(nombre_PD,    "%sPD_lambda%.2f.txt", output_dir, lambda);
        sprintf(nombre_norma, "%snorma_lambda%.2f.txt", output_dir, lambda);
        sprintf(nombre_xesp,  "%sxesp_lambda%.2f.txt", output_dir, lambda);
        sprintf(nombre_Tesp,  "%sTesp_lambda%.2f.txt", output_dir, lambda);
        sprintf(nombre_phiV,  "%sphiV_lambda%.2f.txt", output_dir, lambda);
        sprintf(nombre_PDnd,  "%sPD_nD_todos_lambda%.2f.txt", output_dir, lambda);

        char nombre_Klambda[256];
        sprintf(nombre_Klambda, "%sKvsLambda.txt", output_dir);

        // --- Precalculo matrices fijas ---
        inicializarPotencial(V, lambda, k_tilde);
        construirTridiagonal(A1, A2, A3, V, s_tilde);
        alfaGamma(A1, A2, A3, alfa, gamma);

        // --- Simulaciones Monte Carlo ---
        for (int experimento = 0; experimento < M; experimento++) {
            inicializarPhi(phi, nciclos);

            FILE *f_PD = NULL, *f_norma = NULL, *f_xesp = NULL, *f_Tesp = NULL, *f_phiV = NULL;
            if (experimento == 0) {
                f_PD    = fopen(nombre_PD, "w");
                f_norma = fopen(nombre_norma, "w");
                f_xesp  = fopen(nombre_xesp, "w");
                f_Tesp  = fopen(nombre_Tesp, "w");
                f_phiV  = fopen(nombre_phiV, "w");
                if (!f_PD || !f_norma || !f_xesp || !f_Tesp || !f_phiV) {
                    printf("Error abriendo archivos de guardado\n");
                    return 1;
                }
            }

            for (int n = 0; n < T; n++) {
                PD[n]     = probabilidad(phi);
                normas[n] = norm(phi);
                xesp[n]   = valor_esperado_x(phi);
                Tesp[n]   = valor_esperado_T(phi);

                if (experimento == 0) {
                    fprintf(f_PD,    "%d\t%.10f\n", n, PD[n]);
                    fprintf(f_norma, "%d\t%.10f\n", n, normas[n]);
                    fprintf(f_xesp,  "%d\t%.10f\n", n, xesp[n]);
                    fprintf(f_Tesp,  "%d\t%.10f\n", n, Tesp[n]);
                    if (n % 50 == 0) {
                        for (int j = 0; j <= N; j++) {
                            double x = j * h;
                            double densidad = pow(cabs(phi[j]), 2);
                            fprintf(f_phiV, "%.8f\t%.10f\t%.10f\n", x, densidad, V[j]);
                        }
                        fprintf(f_phiV, "\n");
                    }
                }
                pasoTemporal(phi, chi, A3, b, gamma, alfa, beta, s_tilde, phi_next);
            }

            if (experimento == 0) {
                fclose(f_PD); fclose(f_norma); fclose(f_xesp); fclose(f_Tesp); fclose(f_phiV);
            }

            int nD = maximo(PD, T);
            double PD_nD = PD[nD];
            PD_total += PD_nD;

            FILE *f_PDnd = fopen(nombre_PDnd, "a");
            if (f_PDnd) {
                fprintf(f_PDnd, "%.10f\n", PD_nD);
                fclose(f_PDnd);
            }

            double p = (double)rand() / RAND_MAX;
            if (p < PD_nD) mT++;
        }

        double K = (double)mT / M;
        double PD_media = PD_total / M;
        printf("lambda = %.3f, K = %.10f, PD_media = %.10f\n", lambda, K, PD_media);

        FILE *f_Klambda = fopen(nombre_Klambda, "a");
        if (f_Klambda) {
            fprintf(f_Klambda, "%.3f\t%.10f\t%.10f\n", lambda, K, PD_media);
            fclose(f_Klambda);
        }
    }
    return 0;
}