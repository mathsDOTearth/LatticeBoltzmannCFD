#ifndef LBM_H
#define LBM_H

#define Q 9  // Number of velocity directions in D2Q9

typedef struct {
    int nx;
    int ny;
    double w[Q];
    int cx[Q];
    int cy[Q];
    double ***f;
    double ***feq;
    double **rho;
    double **ux;
    double **uy;
} LBM;

void initialise(LBM *lbm, int nx, int ny);
void compute_macroscopic_variables(LBM *lbm);
void collision_step(LBM *lbm);
void streaming_step(LBM *lbm);
void boundary_conditions(LBM *lbm);
void free_lbm(LBM *lbm);

#endif // LBM_H
