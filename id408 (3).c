#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// ========================== Constants ==========================
#define N 1000000
#define A 16807
#define M 2147483647
#define C 0
#define PI 3.141592653589793

#define Z_MIN 0.0f
#define Z_MAX 200.0f
#define ALBEDO 1.0f

// ========================== Declerations ==========================
long long msg(long long seed);
float invCDF(float y);
float P(float mu);
void rejection();
void cum();
float ray_scatter_step(float *dir, float *r, float phi, float l, long long *seed);
float iso_scatter_step(float *r, float phi, float l, long long *seed);
void scattering(float tau, const char *method);

// ========================== Main ==========================
int main() {
    // Rayleigh Sampling Tests
    rejection();
    cum();  
    // Monte Carlo Slab Scattering
    scattering(10.0f, "iso");
    scattering(10.0f, "ray");
    scattering(0.1f, "ray");

    return 0;
}

// ========================== DEFINITIONS ==========================

// ========================== RNG ==========================
long long msg(long long seed) {return (seed * A + C) % M;}

// ========================== Phase Function Helpers ==========================
float invCDF(float y) {return 1 / cbrt(sqrt(4 * y * y + 1) - 2 * y) - cbrt(sqrt(4 * y * y + 1) - 2 * y);} // Inverse of Cumulative Density Function

float P(float mu) {return 0.375f * (1.0f + mu * mu);  // Rayleigh phase function
}
// ========================== Rejection Sampling ==========================
void rejection() {
    FILE *fp = fopen("Rejection.txt", "w");
    if (!fp) {
        printf("File error!\n");
        return;
    }

    long long seed = 1;
    long double mu = -1.0 + (long double)1.0 / N;
    long double dmu = 2.0 / N;
    int accepted = 0;
    clock_t start = clock();

    for (int i = 0; i < N; i++) {
        seed = msg(seed);
        float y = ((float)seed / M) * 0.75f;

        if (y < P(mu)) {
            fprintf(fp, "%.6Lf %.6f\n", mu, y);
            accepted++;
        }
        mu += dmu;
    }

    clock_t end = clock();
    printf("\nRejection Sampling Efficiency: %.2f%%\n", 100.0f * accepted / N);
    printf("Time: %.3f s\n", (double)(end - start) / CLOCKS_PER_SEC);
    fclose(fp);
}

// ========================== Cumulative Method ==========================
void cum() {
    FILE *fp2 = fopen("Cumulative.txt", "w");
    if (!fp2) {
        printf("File error!\n");
        return;
    }

    long long seed = 1;
    int accepted = 0;
    int N_bins = 50;
    float dmu = 2.0f / N_bins;
    int binned[N_bins];
    for (int i = 0; i < N_bins; i++) {
        binned[i] = 0;
    }

    clock_t start = clock();

    for (int i = 0; i < N; i++) {
        seed = msg(seed);
        double y = (double)seed / M;

        float mu = invCDF(y);
        seed = msg(seed);
        if (seed % 2 == 0) mu = -mu;

        int index = (int)((mu + 1.0f) * N_bins / 2);
            binned[index]++;
            accepted++;
    }

    clock_t end = clock();

    for (int i = 0; i < N_bins; i++) {
        float mu_center = -1.0f + (i + 0.5f) * dmu;
        float density = binned[i] / (float)(N * dmu);
        fprintf(fp2, "%.6f %.6f %.6f \n", mu_center, density,dmu);
    }

    printf("\nInverse Transform Sampling Efficiency: %.2f%%\n", 100.0f * accepted / N);
    printf("Time: %.3f s\n", (double)(end - start) / CLOCKS_PER_SEC);
    fclose(fp2);
}
// ========================== Scattering Functions ==========================
float ray_scatter_step(float *dir, float *r, float phi, float l, long long *seed) {
    *seed = msg(*seed);
    float u = (float)(*seed) / M;
    float mu = invCDF(u);

    *seed = msg(*seed);
    if (*seed % 2 == 0)
        mu = -mu;

    float theta = acosf(mu);

    float perp[3];
    if (fabsf(dir[2]) > 0.999f) {
        perp[0] = 1; perp[1] = 0; perp[2] = 0;
    } else {
        perp[0] = -dir[1];
        perp[1] = dir[0];
        perp[2] = 0;
    }

    float norm = sqrtf(perp[0]*perp[0] + perp[1]*perp[1] + perp[2]*perp[2]);
    float n[3] = {perp[0]/norm, perp[1]/norm, perp[2]/norm};

    float v[3] = {
        dir[1]*n[2] - dir[2]*n[1],
        dir[2]*n[0] - dir[0]*n[2],
        dir[0]*n[1] - dir[1]*n[0]
    };

    float d_new[3] = {
        sinf(theta)*cosf(phi)*n[0] + sinf(theta)*sinf(phi)*v[0] + cosf(theta)*dir[0],
        sinf(theta)*cosf(phi)*n[1] + sinf(theta)*sinf(phi)*v[1] + cosf(theta)*dir[1],
        sinf(theta)*cosf(phi)*n[2] + sinf(theta)*sinf(phi)*v[2] + cosf(theta)*dir[2]
    };

    norm = sqrtf(d_new[0]*d_new[0] + d_new[1]*d_new[1] + d_new[2]*d_new[2]);
    for (int j = 0; j < 3; j++) {
        dir[j] = d_new[j] / norm;
        r[j] += l * dir[j];
    }

    return theta;
}

float iso_scatter_step(float *r, float phi, float l, long long *seed) {
    *seed = msg(*seed);
    float u = (float)(*seed) / M;

    float theta = acosf(1.0f - 2.0f * u);

    float dx = l * sinf(theta) * cosf(phi);
    float dy = l * sinf(theta) * sinf(phi);
    float dz = l * cosf(theta);

    r[0] += dx;
    r[1] += dy;
    r[2] += dz;
    return theta;
}

// ========================== Monte Carlo Slab Scattering ==========================
void scattering(float tau, const char *method) {
    const float ALPHA_NU = tau / (Z_MAX - Z_MIN);

    char filename[128];
    snprintf(filename, sizeof(filename), "tau_%.1f_%s.txt", tau, method);

    FILE *f = fopen(filename, "w");
    if (!f) {
        printf("Error opening file '%s'.\n", filename);
        return;
    }

    long long seed = 1;
    int absorbed = 0;
    int N_bins = 10;
    int binned[N_bins];
    for (int i = 0; i < N_bins; i++) {
        binned[i] = 0;
    }
    float dmu=N_bins / 2.0f;
    for (int i = 0; i < N; i++) {
        float r[3] = {0.0f, 0.0f, 0.0f};
        float dir[3] = {0.0f, 0.0f, 1.0f};

        do {
            seed = msg(seed);
            float zeta = (float)seed / M;

            if (zeta < ALBEDO) {
                seed = msg(seed);
                float y = (float)seed / M;
                float tau_sample = -logf(y);
                float l = tau_sample / ALPHA_NU;

                seed = msg(seed);
                float phi = 2.0f * PI * ((float)seed / M);

                float theta;
                if (strcmp(method, "ray") == 0)
                    theta = ray_scatter_step(dir, r, phi, l, &seed);
                else
                    theta = iso_scatter_step(r, phi, l, &seed);

                if (r[2] > Z_MAX) {
                    float mu_out = cosf(theta);
                    int index = (int)((mu_out + 1.0f) * dmu);
                    if (index >= 0 && index < N_bins)
                        binned[index]++;
                    break;
                }
            } else {
                absorbed++;
                break;
            }
        } while (r[2] >= Z_MIN && r[2] <= Z_MAX);
    }

    for (int i = 0; i < N_bins; i++) {
        float mu = -1.0f + (i + 0.5f) * (dmu);
        float intensity = (float)binned[i] / N;
        fprintf(f, "%f\t%f %f\n", mu, intensity,dmu);
    }

    fclose(f);
    printf("Simulated tau=%.1f, method=%s → %s | Absorbed: %d (%.2f%%)\n",
           tau, method, filename, absorbed, 100.0f * absorbed / N);
}
