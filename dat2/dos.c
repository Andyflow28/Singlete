#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Qromb.h"
#include "Utils.h"

#define N 1000
#define DELTA0 1.0 // Define Delta0 como una constante global
#define M_PI 3.14159265358979323846

float A[N], B[N];
double En = -0.4;
double Hop_t = -0.2;
int current_index; // Variable global para almacenar el índice actual

// Declaraciones de funciones
double Func1(double x);
double Func2(double x);
double dos_f(double x);

int main(void) {
    double LowerLimitDx = M_PI - acos(En / (4 * Hop_t));
    double normalizacion, f, g, DOS;

    char* file_names[] = {"g0001delta0=1.dat", "g0005delta0=1.dat", "g0010delta0=1.dat", "g0015delta0=1.dat", "g0020delta0=1.dat"};
    char* DOS_names[] = {"DOSc0g0001.dat", "DOSc0g0005.dat", "DOSc0g0010.dat", "DOSc0g0015.dat", "DOSc0g0020.dat"};

    for (int i = 0; i < 5; i++) {
        FILE *fp = fopen(file_names[i], "r");
        if (!fp) {
            perror("Error al abrir archivo de entrada");
            exit(EXIT_FAILURE);
        }
        for (int j = 0; j < N; j++) {
            fscanf(fp, "%f %f", &A[j], &B[j]);
        }
        fclose(fp);

        FILE *dos = fopen(DOS_names[i], "w");
        if (!dos) {
            perror("Error al abrir archivo de salida");
            exit(EXIT_FAILURE);
        }

        for (int n = 0; n < N; n++) {
            current_index = n; // Actualizamos el índice global
            normalizacion = qromb(dos_f, LowerLimitDx, M_PI);
            f = qromb(Func1, LowerLimitDx, M_PI);
            g = qromb(Func2, LowerLimitDx, M_PI);
            DOS = (f + g) / normalizacion;
            fprintf(dos, "%f\t%f\n", A[n], DOS);
        }
        fclose(dos);
    }
    return 0;
}

double Func1(double x) {
    int n = current_index; // Usamos el índice global
    double cos_x = cos(x);
    double sin_x = sqrt(1 - cos_x * cos_x);
    double aux1 = En + 2 * Hop_t * cos_x;
    double z = M_PI - acos(0.5 * aux1 / Hop_t);
    double cos_z = cos(z);
    double sin_z = sqrt(1 - cos_z * cos_z);
    double vx = 2 * Hop_t * sin_x;
    double vz = 2 * Hop_t * sin_z;
    double vel = sqrt(vx * vx + vz * vz);
    double jacobian = sqrt(1 + pow(vx / vz, 2)) / vel;

    double a_k = pow(A[n], 2.0) - pow(B[n], 2.0) - pow(cos_x - cos_z, 2.0);
    double b = 2 * A[n] * B[n];
    double rho_k = sqrt(pow(a_k, 2.0) + pow(b, 2.0));
    double positive_k = 1 + a_k / rho_k;
    double result = jacobian * sqrt(2 * rho_k) * (A[n] / sqrt(2 * rho_k)) * sqrt(positive_k);
    return result;
}

double Func2(double x) {
    int n = current_index; // Usamos el índice global
    // Implementar lógica similar a Func1 si es necesario
    return 0.0;
}

double dos_f(double x) {
    double cos_x = cos(x);
    double sin_x = sqrt(1 - cos_x * cos_x);
    double aux1 = En + 2 * Hop_t * cos_x;
    double z = M_PI - acos(0.5 * aux1 / Hop_t);
    double vx = 2 * Hop_t * sin_x;
    double vz = 2 * Hop_t * sin(z);
    double vel = sqrt(vx * vx + vz * vz);
    double jacobian = sqrt(1 + pow(vx / vz, 2)) / vel;
    return jacobian;
}
