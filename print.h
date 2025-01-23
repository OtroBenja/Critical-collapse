#include <time.h>
#include <stdio.h>
#include "constants.h"

void print_mass(double *MassHist,struct tm time_data,int print_iterations){
    char fileName[50];
    snprintf(fileName, sizeof(fileName), "Mass_%02d%02d%02d.dat", time_data.tm_hour, time_data.tm_min, time_data.tm_sec);
    FILE* data = fopen(fileName,"w");
    //Print M
    for(int i=0;i<print_iterations-1;i++){
        fprintf(data,"%lf,",MassHist[i]);
    }
    fprintf(data,"%lf",MassHist[print_iterations-1]);
    fclose(data);
}

void print_data(float **hist,int fType,double *model_parameters,int iterations,double maxR,double deltaR,int nGrids,time_t totalTime){
    int print_iterations, printR;
    if(SAVE_MODE == 0){
        print_iterations = iterations/SAVE_ITERATION;
        printR = (int)((maxR/deltaR)/SAVE_RES);
    }
    if(SAVE_MODE == 1){
        print_iterations = iterations;
        printR = (int)((fmin(maxR,MAX_R)-MIN_R)/deltaR);
    }
    double p0 = model_parameters[0];
    double r0 = model_parameters[1];
    double d  = model_parameters[2];
    //Add time to filename
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    //if(MASS) print_mass(hist[6],tm,print_iterations);
    char fileName[50];
    snprintf(fileName, sizeof(fileName), "Output_%02d%02d%02d.dat", tm.tm_hour, tm.tm_min, tm.tm_sec);
    FILE* data = fopen(fileName,"w");
    //Print all parameters

    if(nGrids == 1){
             if(METRIC == 0) fprintf(data,"Metric: Minkowski\n");
        else if(METRIC == 1) fprintf(data,"Metric: Choptuik\n");
        else if(METRIC == 2) fprintf(data,"Metric: Quasi-Static Choptuik\n");
        else                 fprintf(data,"Metric: Other\n");
    }
    if(nGrids > 1){
             if(SUBGRID_MODE == 0) fprintf(data,"Metric: Variable Choptuik A\n");
        else if(SUBGRID_MODE == 1) fprintf(data,"Metric: Variable Choptuik B\n");
    }
         if(fType == 0) fprintf(data,"Function type: Hyperbolic Tan\n");
    else if(fType == 1) fprintf(data,"Function type: Exponential\n");
    else if(fType == 2) fprintf(data,"Function type: Constant\n");
    else if(fType == 3) fprintf(data,"Function type: Linear\n");
    fprintf(data,"p0: %.10e\n",p0);
    fprintf(data,"r0: %lf\n",r0);
    fprintf(data,"d: %lf\n",d);
    fprintf(data,"R step size: %lf\n",deltaR);
    fprintf(data,"Maximum R: %lf\n",maxR);
    fprintf(data,"Final BH radius: %e\n",bh_radius);
    fprintf(data,"Final BH mass: %e\n",bh_mass);
    fprintf(data,"Iterations: %d\n",iterations);
    fprintf(data,"Number of grids: %d\n",nGrids);
    fprintf(data,"Total simulation time: %ld\n",totalTime);
    //Print R
    for(int ir=0;ir<(printR-1);ir++){
        fprintf(data,"%lf,",hist[0][ir]);
    }
    fprintf(data,"%lf\n",hist[0][printR-1]);
    //Print phi
    for(int i=0;i<print_iterations;i++){
        for(int ir=0;ir<(printR-1);ir++){
            fprintf(data,"%e,",hist[1][i*printR+ir]);
        }
        fprintf(data,"%e\n",hist[1][i*printR+printR-1]);
    }
    //Print mass
    for(int i=0;i<print_iterations;i++){
        for(int ir=0;ir<(printR-1);ir++){
            fprintf(data,"%lf,",hist[6][i*printR+ir]);
        }
        fprintf(data,"%lf\n",hist[6][i*printR+printR-1]);
    }
    //Print Phi
    for(int i=0;i<print_iterations;i++){
        for(int ir=0;ir<(printR-1);ir++){
            fprintf(data,"%e,",hist[2][i*printR+ir]);
        }
        fprintf(data,"%e\n",hist[2][i*printR+printR-1]);
    }
    //Print Pi
    for(int i=0;i<print_iterations;i++){
        for(int ir=0;ir<(printR-1);ir++){
            fprintf(data,"%e,",hist[3][i*printR+ir]);
        }
        fprintf(data,"%e\n",hist[3][i*printR+printR-1]);
    }
    //Print a
    for(int i=0;i<print_iterations;i++){
        for(int ir=0;ir<(printR-1);ir++){
            fprintf(data,"%lf,",hist[4][i*printR+ir]);
        }
        fprintf(data,"%lf\n",hist[4][i*printR+printR-1]);
    }
    //Print alpha
    for(int i=0;i<print_iterations-1;i++){
        for(int ir=0;ir<(printR-1);ir++){
            fprintf(data,"%lf,",hist[5][i*printR+ir]);
        }
        fprintf(data,"%lf\n",hist[5][i*printR+printR-1]);
    }
    for(int ir=0;ir<(printR-1);ir++){
        fprintf(data,"%lf,",hist[5][(print_iterations-1)*printR+ir]);
    }
    fprintf(data,"%lf",hist[5][print_iterations*printR-1]);
    fclose(data);
}