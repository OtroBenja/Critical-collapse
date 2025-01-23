#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#if defined(_OPENMP)
#include <omp.h>
#else 
int omp_get_num_threads(){return 1;}
int omp_get_max_threads(){return 1;}
int omp_get_thread_num(){return 1;}
void omp_set_num_threads(){return;}
#endif

double max_a;
double bh_mass;
double bh_radius;
bool is_bh;
int last_iteration;
double maxR_global;
double norm_global;

#include "initialize.h"
#include "print.h"
#include "variable_iteration.h"
#include "variable_integration.h"



void check_collapse(double ***all_subgrids, int grid_n, int iteration){

    double **grid_l;
    double *a;
    int nR;

    for(int ig=N_LEVELS;ig>grid_n-1;ig--){
        //Load each subgrid (in order)
        grid_l = all_subgrids[ig];
        a = grid_l[4];
        nR = (int)(*(grid_l[8]));

        //Check if a collapse has happened if a > 1/TOLERANCE^0.5
        for(int ir=0;ir<nR;ir++){
            if(a[ir]>max_a){
                //If a is over the threshold, find the local maxima
                for(;ir<nR;){
                    if(a[ir+1]>a[ir]) ir++;
                    else break;
                }
                double *r = grid_l[0];
                bh_radius = r[ir];
                bh_mass = mass_ir(all_subgrids,ig,ir);
                last_iteration = iteration;
                is_bh = true;
                break;
            }
        }
    }
}


double **initialize_subgrid(double fType,double *model_params,double initial_r, double final_r, double deltaR){
    double **initial_field;

    // If initial and final r are not in the boundaries of the simulation, add ghost zones

    double new_initial_r;
    int  left_ghost_size;
    bool left_limit = false;
    if (initial_r == 0) {
        new_initial_r = initial_r;
        left_ghost_size = 0;
        left_limit = true;
    } else {
        new_initial_r = initial_r - GHOST_SIZE*deltaR;
        left_ghost_size = GHOST_SIZE;
    }

    double   new_final_r;
    int right_ghost_size;
    bool right_limit = false;
    if (final_r == maxR_global) {
        new_final_r = final_r + deltaR*0.01;
        right_ghost_size = 0;
        right_limit = true;
    } else {
        new_final_r = final_r + GHOST_SIZE*deltaR + deltaR*0.01;
        right_ghost_size = GHOST_SIZE;
    }

    initial_field = initialize_field(fType,model_params,deltaR,new_initial_r,new_final_r);

    int nR = (int)((deltaR*0.01+final_r-initial_r)/deltaR);
    double nR_f = (double)nR;
    int nR_ext = (int)((new_final_r-new_initial_r)/deltaR);
    //int nR_ext = nR + left_ghost_size + right_ghost_size;
    double deltaT = deltaR/5.;
    double *a     = malloc(sizeof(double)*(nR_ext));
    double *alpha = malloc(sizeof(double)*(nR_ext));
    for(int ir=0; ir<(nR_ext); ir++){
        a[ir] = 1.0;
        alpha[ir] = 1.0;
    }
    double *nR_p = malloc(sizeof(double));
    double *dR_p = malloc(sizeof(double));
    double *dT_p = malloc(sizeof(double));
    double *ri_p = malloc(sizeof(double));
    double *rf_p = malloc(sizeof(double));
    double *llim_p = malloc(sizeof(double));
    double *rlim_p = malloc(sizeof(double));
    double *partial_a_p = malloc(sizeof(double));
    double *partial_alpha_p = malloc(sizeof(double));

    dR_p[0] = deltaR;
    dT_p[0] = deltaT;
    nR_p[0] = nR_f;
    ri_p[0] = initial_r;
    llim_p[0] = (double)left_limit;
    rlim_p[0] = (double)right_limit;
    partial_a_p[0] = nR_f;
    partial_alpha_p[0] = nR_f;

    if(left_limit) 
        initial_field[0][left_ghost_size] = 1.0E-50;
    rf_p[0] = final_r-deltaR;

    printf("nR: %d\n",nR);
    printf("nR_extended: %d\n",nR_ext);
    printf("dR_pointer: %lf\n",*dR_p);
    printf("dT_pointer: %lf\n",*dT_p);
    //printf("ri: %.20lf\n",initial_r);
    //printf("rf: %.20lf\n",final_r);
    printf("grid range: %e - %e\n\n",initial_field[0][left_ghost_size],(initial_field[0]+left_ghost_size)[nR-1]);

    double **subgrid = malloc(sizeof(double*)*15);
    subgrid[0]  = initial_field[0]+left_ghost_size; // r values of the grid
    subgrid[1]  = initial_field[1]+left_ghost_size; // phi values of the grid
    subgrid[2]  = initial_field[2]+left_ghost_size; // Phi values of the grid
    subgrid[3]  = initial_field[3]+left_ghost_size; // Pi values of the grid
    subgrid[4]  =     a+left_ghost_size; // a values of the grid
    subgrid[5]  = alpha+left_ghost_size; // alpha values of the grid
    subgrid[6]  = dR_p;   // deltaR of the grid (as pointer)
    subgrid[7]  = dT_p;   // deltaT of the grid (as pointer)
    subgrid[8]  = nR_p;   // number of points in the grid (as pointer)
    subgrid[9]  = ri_p;   // min r of the grid (as pointer)
    subgrid[10] = rf_p;   // max r of the grid (as pointer)
    subgrid[11] = llim_p; // bool for left limit of the grid (as pointer)
    subgrid[12] = rlim_p; // bool for right limit of the grid (as pointer)
    subgrid[13] = partial_a_p;     //Partial sum of a in the subgrid (pre-normalization)
    subgrid[14] = partial_alpha_p; //Partial sum of alpha in the subgrid (pre-normalization)

    return subgrid;
}


int main(int argc, char* argv[]){
    double* model_parameters = malloc(sizeof(double)*3);
    double** initial_conditions;
    float** hist;
    double* r;
    double* phi;
    double* Phi;
    double* Pi;
    //Define simulation parameters
    int fType = 0;
    if((argc>1) && atoi(argv[1])) fType = atoi(argv[1]);
    double p0 = 0.000008;
    if((argc>2) && atof(argv[2])) p0 = atof(argv[2]);
    double r0 = 20.;
    if((argc>3) && atof(argv[3])) r0 = atof(argv[3]);
    double d = 3.;
    if((argc>4) && atof(argv[4])) d  = atof(argv[4]);
    model_parameters[0] = p0;
    model_parameters[1] = r0;
    model_parameters[2] = d;
    
    //Define simulation limits
    double deltaR = 0.01;
    if((argc>5) && atof(argv[5])) deltaR = atof(argv[5]);
    double maxR = 50;
    if((argc>6) && atoi(argv[6])) maxR = atof(argv[6]);
    //maxR_global = maxR;
    maxR_global = maxR+deltaR*((double)EXTRA_N);
    int iterations = ITERATIONS;
    if((argc>7) && atoi(argv[7])) iterations = atoi(argv[7]);
    printf("Total iterations: %d\n",iterations);

    //Pass initial conditions to iteration
    time_t initTime = time(NULL);
    hist = variable_iteration(fType, model_parameters,deltaR,maxR,iterations,SAVE_ITERATION);
    time_t finalTime = time(NULL);
    time_t timeDelta = (finalTime-initTime);

    //Print simulation history to a file
    if(SAVE_MODE == 0)
        print_data(hist,fType,model_parameters,last_iteration,maxR,deltaR,N_LEVELS,timeDelta);
    else if(SAVE_MODE == 1){
        print_data(hist,fType,model_parameters,last_iteration-FIRST_ITERATION,maxR,deltaR,N_LEVELS,timeDelta);
        }
    printf("Finished, total time: %lds\n", timeDelta);
}

