#include "lbm.h"
#include <stdlib.h>
#include <stdlib.h>  // Include for rand()
#include <time.h>    // Include for time()

void initialise(LBM *lbm, int nx, int ny) {
    lbm->nx = nx;
    lbm->ny = ny;

    srand(time(NULL)); // Seed the random number generator

    double w[Q] = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
    int cx[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    int cy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

    for (int i = 0; i < Q; i++) {
        lbm->w[i] = w[i];
        lbm->cx[i] = cx[i];
        lbm->cy[i] = cy[i];
    }

    lbm->f = malloc(nx * sizeof(double**));
    lbm->feq = malloc(nx * sizeof(double**));
    lbm->rho = malloc(nx * sizeof(double*));
    lbm->ux = malloc(nx * sizeof(double*));
    lbm->uy = malloc(nx * sizeof(double*));

    for (int x = 0; x < nx; x++) {
        lbm->f[x] = malloc(ny * sizeof(double*));
        lbm->feq[x] = malloc(ny * sizeof(double*));
        lbm->rho[x] = malloc(ny * sizeof(double));
        lbm->ux[x] = malloc(ny * sizeof(double));
        lbm->uy[x] = malloc(ny * sizeof(double));

        for (int y = 0; y < ny; y++) {
            lbm->f[x][y] = malloc(Q * sizeof(double));
            lbm->feq[x][y] = malloc(Q * sizeof(double));

            lbm->rho[x][y] = 1.0;

            // Increase the magnitude of the initial velocity perturbation
            lbm->ux[x][y] = 0.5 * ((double)rand() / RAND_MAX - 0.5);
            lbm->uy[x][y] = 0.5 * ((double)rand() / RAND_MAX - 0.5);

            for (int k = 0; k < Q; k++) {
                lbm->feq[x][y][k] = lbm->w[k] * lbm->rho[x][y] * (1 + 3 * (lbm->cx[k] * lbm->ux[x][y] + lbm->cy[k] * lbm->uy[x][y]));
                lbm->f[x][y][k] = lbm->feq[x][y][k];
            }
        }
    }
}

void compute_macroscopic_variables(LBM *lbm) {
    for (int x = 0; x < lbm->nx; x++) {
        for (int y = 0; y < lbm->ny; y++) {
            double sum_f = 0.0;
            double sum_ux = 0.0;
            double sum_uy = 0.0;

            for (int k = 0; k < Q; k++) {
                sum_f += lbm->f[x][y][k];
                sum_ux += lbm->f[x][y][k] * lbm->cx[k];
                sum_uy += lbm->f[x][y][k] * lbm->cy[k];
            }

            lbm->rho[x][y] = sum_f;
            lbm->ux[x][y] = sum_ux / lbm->rho[x][y];
            lbm->uy[x][y] = sum_uy / lbm->rho[x][y];
        }
    }
}

void collision_step(LBM *lbm) {
    double tau = 0.6;  // Relaxation time

    for (int x = 0; x < lbm->nx; x++) {
        for (int y = 0; y < lbm->ny; y++) {
            for (int k = 0; k < Q; k++) {
                lbm->feq[x][y][k] = lbm->w[k] * lbm->rho[x][y] * (1 + 3 * (lbm->cx[k] * lbm->ux[x][y] + lbm->cy[k] * lbm->uy[x][y]));
                lbm->f[x][y][k] = lbm->f[x][y][k] - (lbm->f[x][y][k] - lbm->feq[x][y][k]) / tau;
            }
        }
    }
}

void streaming_step(LBM *lbm) {
    // Allocate memory on the heap for the temporary array
    double ***temp = malloc(lbm->nx * sizeof(double**));
    for (int x = 0; x < lbm->nx; x++) {
        temp[x] = malloc(lbm->ny * sizeof(double*));
        for (int y = 0; y < lbm->ny; y++) {
            temp[x][y] = malloc(Q * sizeof(double));
        }
    }

    // Copy the current state to the temporary array
    for (int x = 0; x < lbm->nx; x++) {
        for (int y = 0; y < lbm->ny; y++) {
            for (int k = 0; k < Q; k++) {
                temp[x][y][k] = lbm->f[x][y][k];
            }
        }
    }

    // Perform streaming
    for (int x = 0; x < lbm->nx; x++) {
        for (int y = 0; y < lbm->ny; y++) {
            for (int k = 0; k < Q; k++) {
                int x_new = (x + lbm->cx[k] + lbm->nx) % lbm->nx;
                int y_new = (y + lbm->cy[k] + lbm->ny) % lbm->ny;
                lbm->f[x_new][y_new][k] = temp[x][y][k];
            }
        }
    }

    // Free the allocated memory
    for (int x = 0; x < lbm->nx; x++) {
        for (int y = 0; y < lbm->ny; y++) {
            free(temp[x][y]);
        }
        free(temp[x]);
    }
    free(temp);
}

void boundary_conditions(LBM *lbm) {
    // Implement boundary conditions (e.g., bounce-back for solid walls)
}

void free_lbm(LBM *lbm) {
    for (int x = 0; x < lbm->nx; x++) {
        for (int y = 0; y < lbm->ny; y++) {
            free(lbm->f[x][y]);
            free(lbm->feq[x][y]);
        }
        free(lbm->f[x]);
        free(lbm->feq[x]);
        free(lbm->rho[x]);
        free(lbm->ux[x]);
        free(lbm->uy[x]);
    }
    free(lbm->f);
    free(lbm->feq);
    free(lbm->rho);
    free(lbm->ux);
    free(lbm->uy);
}
