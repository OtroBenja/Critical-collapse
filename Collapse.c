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

double bh_mass;
double bh_radius;
bool is_bh;
int last_iteration;

#include "initialize.h"
#include "iteration.h"
#include "print.h"

int main(int argc, char* argv[]){
    double* model_parameters = malloc(sizeof(double)*3);
    double** initial_conditions;
    double** hist;
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

    //Set Kreiss Oliger coefficient for dampening
    double epsilon = EPSILON;
    if((argc>5) && atof(argv[5])) epsilon  = atof(argv[5]);
    
    //Define simulation limits
    double deltaR = 0.01;
    if((argc>6) && atof(argv[6])) deltaR = atof(argv[6]);
    double maxR = 50;
    if((argc>7) && atoi(argv[7])) maxR = atof(argv[7]);
    int iterations = ITERATIONS;
    if((argc>8) && atoi(argv[8])) iterations = atoi(argv[8]);
    printf("Total iterations: %d\n",iterations);

    initial_conditions = initialize_field(fType,model_parameters,deltaR,0.0,maxR);
    r = initial_conditions[0];
    phi = initial_conditions[1];
    Phi = initial_conditions[2];
    Pi = initial_conditions[3];
    //Pass initial conditions to iteration
    omp_set_num_threads(1);
    if((argc>9) && atoi(argv[9])) omp_set_num_threads(atoi(argv[9]));
    time_t initTime = time(NULL);
    hist = iteration(r,phi,Phi,Pi,deltaR,maxR,iterations,SAVE_ITERATION);
    time_t finalTime = time(NULL);
    int nP = omp_get_max_threads();
    time_t timeDelta = (finalTime-initTime);

    //Print simulation history to a file
    if(SAVE_MODE == 0)
        print_data(hist,fType,model_parameters,last_iteration,maxR,deltaR,nP,timeDelta,epsilon);
    else if(SAVE_MODE == 1){
        print_data(hist,fType,model_parameters,last_iteration-FIRST_ITERATION,maxR-MIN_R,deltaR,nP,timeDelta,epsilon);
        }
    printf("Finished, total time: %lds\n", timeDelta);
}




