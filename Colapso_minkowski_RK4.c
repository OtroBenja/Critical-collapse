#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SAVE_RES 100
#define SAVE_ITERATION 100
#define ITERATIONS 100000
#define PI 3.1415



double** initialize_field(double p0,double r0,double d,double q,double deltaR,int maxR){
    int nR = maxR/deltaR;
    double* r = malloc(sizeof(double)*nR);
    double* phi = malloc(sizeof(double)*nR);
    double* Phi = malloc(sizeof(double)*nR);
    double* Pi = malloc(sizeof(double)*nR);
    double** return_values = malloc(sizeof(double*)*4);

    //Define r and calculate initial phi
    for(int i=0;i<nR;i++){
        r[i] = i*deltaR;
        phi[i] = p0*tanh((i*deltaR-r0)/d);
    }

    //Calculate initial Phi
    Phi[0] = 0;
    Phi[nR-1] = 0;
    for(int i=1;i<nR-1;i++){
        Phi[i] = 0.5*(phi[i+1]-phi[i-1])/deltaR;
    }

    //Set initial Pi to zero
    for(int i=0;i<nR;i++){
        Pi[i] = 0;
    }
        
    return_values[0] = r;
    return_values[1] = phi;
    return_values[2] = Phi;
    return_values[3] = Pi;
 
    return return_values;
}

double** iteration(double* r,double* phi,double* Phi,double* Pi,double deltaR,int maxR,int iterations,int save_iteration){

    int nR = maxR/deltaR;
    double deltaT=deltaR/5.;
    double* Pi_t = malloc(sizeof(double)*nR);
    double* Phi_t = malloc(sizeof(double)*nR);
    double* temp1 = malloc(sizeof(double)*nR);
    double* temp2 = malloc(sizeof(double)*nR);
    double* k1 = malloc(sizeof(double)*nR);
    double* l1 = malloc(sizeof(double)*nR);
    double* k2 = malloc(sizeof(double)*nR);
    double* l2 = malloc(sizeof(double)*nR);
    double* k3 = malloc(sizeof(double)*nR);
    double* l3 = malloc(sizeof(double)*nR);
    double* k4 = malloc(sizeof(double)*nR);
    double* l4 = malloc(sizeof(double)*nR);
    double* Xhistory = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/save_iteration));
    double* Yhistory = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/save_iteration));
    double* Fhistory = malloc(sizeof(double)*(nR/SAVE_RES)*(iterations/save_iteration));
    double* Rhistory = malloc(sizeof(double)*(nR/SAVE_RES));
    double** hist = malloc(sizeof(double*)*4);
    int save_count = save_iteration;

    double a = 1.;
    double alpha = 1.;

    for(int i=0;i<iterations;i++){
        //Solve a and alpha with known values of PHI and PI
        //Save values of X and Y
        if(save_count == save_iteration){
            printf("iteration %d\n",i);
            for(int ir=0;ir<(nR/SAVE_RES);ir++){
                //Xhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = r[ir*SAVE_RES]*Phi[ir*SAVE_RES]*sqrt(2*PI)/a;
                //Yhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = r[ir*SAVE_RES]* Pi[ir*SAVE_RES]*sqrt(2*PI)/a;
                Xhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = Phi[ir*SAVE_RES];
                Yhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] =  Pi[ir*SAVE_RES];
                Fhistory[(i/save_iteration)*(nR/SAVE_RES)+(ir)] = phi[ir*SAVE_RES];
            }
            save_count=0;
        }
        save_count+=1;
        //Advance Pi and Phi using RK4
        //calculate k1 and l1
        k1[0] = 0.5*(-3*Pi[0] +4*Pi[1] -Pi[2])/deltaT;
        k1[nR-1] = -0.5*(3*Phi[nR-1] -4*Phi[nR-2] +Phi[nR-3])/deltaT;
        l1[0] = 0.5*(-3*r[0]*r[0]*Phi[0] +4*r[1]*r[1]*Phi[1] -r[2]*r[2]*Phi[2])/(deltaT*(r[0]*r[0]+0.01*deltaR));
        l1[nR-1] = -0.5*(3*r[nR-1]*r[nR-1]*Pi[nR-1] -4*r[nR-2]*r[nR-2]*Pi[nR-2] +r[nR-3]*r[nR-3]*Pi[nR-3])/(deltaT*r[nR-1]*r[nR-1]);
        for(int ir=1;ir<nR-1;ir++){
            k1[ir] = 0.5*(alpha*Pi[ir+1]/a -alpha*Pi[ir-1]/a)/deltaT;
            l1[ir] = 0.5*(r[ir+1]*r[ir+1]*alpha*Phi[ir+1]/a -r[ir-1]*r[ir-1]*alpha*Phi[ir-1]/a)/(deltaT*r[ir]*r[ir]);
        }
        //calculate k2 and l2
        k2[0] = 0.5*(-3*(Pi[0]+0.5*l1[0]*deltaT) +4*(Pi[1]+0.5*l1[1]*deltaT) -(Pi[2]+0.5*l1[2]*deltaT))/deltaT;
        k2[nR-1] = -0.5*(3*(Phi[nR-1]+0.5*k1[nR-1]*deltaT) -4*(Phi[nR-2]+0.5*k1[nR-2]*deltaT) +(Phi[nR-3]+0.5*k1[nR-3]*deltaT))/deltaT;
        l2[0] = 0.5*(-3*r[0]*r[0]*(Phi[0]+0.5*k1[0]*deltaT) +4*r[1]*r[1]*(Phi[1]+0.5*k1[1]*deltaT) -r[2]*r[2]*(Phi[2]+0.5*k1[2]*deltaT))/(deltaT*(r[0]*r[0]+0.01*deltaR));
        l2[nR-1] = -0.5*(3*r[nR-1]*r[nR-1]*(Pi[nR-1]+0.5*l1[nR-1]*deltaT) -4*r[nR-2]*r[nR-2]*(Pi[nR-2]+0.5*l1[nR-2]*deltaT) +r[nR-3]*r[nR-3]*(Pi[nR-3]+0.5*l1[nR-3]*deltaT))/(deltaT*r[nR-1]*r[nR-1]);
        for(int ir=1;ir<nR-1;ir++){
            k2[ir] = 0.5*(alpha*(Pi[ir+1]+0.5*l1[ir+1]*deltaT)/a -alpha*(Pi[ir-1]+0.5*l1[ir-1]*deltaT)/a)/deltaT;
            l2[ir] = 0.5*(r[ir+1]*r[ir+1]*alpha*(Phi[ir+1]+0.5*k1[ir+1]*deltaT)/a -r[ir-1]*r[ir-1]*alpha*(Phi[ir-1]+0.5*k1[ir-1]*deltaT)/a)/(deltaT*r[ir]*r[ir]);
        }
        //calculate k3 and l3
        k2[0] = 0.5*(-3*(Pi[0]+0.5*l2[0]*deltaT) +4*(Pi[1]+0.5*l2[1]*deltaT) -(Pi[2]+0.5*l2[2]*deltaT))/deltaT;
        k2[nR-1] = -0.5*(3*(Phi[nR-1]+0.5*k2[nR-1]*deltaT) -4*(Phi[nR-2]+0.5*k2[nR-2]*deltaT) +(Phi[nR-3]+0.5*k2[nR-3]*deltaT))/deltaT;
        l2[0] = 0.5*(-3*r[0]*r[0]*(Phi[0]+0.5*k2[0]*deltaT) +4*r[1]*r[1]*(Phi[1]+0.5*k2[1]*deltaT) -r[2]*r[2]*(Phi[2]+0.5*k2[2]*deltaT))/(deltaT*(r[0]*r[0]+0.01*deltaR));
        l2[nR-1] = -0.5*(3*r[nR-1]*r[nR-1]*(Pi[nR-1]+0.5*l2[nR-1]*deltaT) -4*r[nR-2]*r[nR-2]*(Pi[nR-2]+0.5*l2[nR-2]*deltaT) +r[nR-3]*r[nR-3]*(Pi[nR-3]+0.5*l2[nR-3]*deltaT))/(deltaT*r[nR-1]*r[nR-1]);
        for(int ir=1;ir<nR-1;ir++){
            k3[ir] = 0.5*(alpha*(Pi[ir+1]+0.5*l2[ir+1]*deltaT)/a -alpha*(Pi[ir-1]+0.5*l2[ir-1]*deltaT)/a)/deltaT;
            l3[ir] = 0.5*(r[ir+1]*r[ir+1]*alpha*(Phi[ir+1]+0.5*k2[ir+1]*deltaT)/a -r[ir-1]*r[ir-1]*alpha*(Phi[ir-1]+0.5*k2[ir-1]*deltaT)/a)/(deltaT*r[ir]*r[ir]);
        }
        //calculate k4 and l4
        k4[0] = 0.5*(-3*(Pi[0]+l3[0]*deltaT) +4*(Pi[1]+l3[1]*deltaT) -(Pi[2]+l3[2]*deltaT))/deltaT;
        k4[nR-1] = -0.5*(3*(Phi[nR-1]+k3[nR-1]*deltaT) -4*(Phi[nR-2]+k3[nR-2]*deltaT) +(Phi[nR-3]+k3[nR-3]*deltaT))/deltaT;
        l4[0] = 0.5*(-3*r[0]*r[0]*(Phi[0]+k3[0]*deltaT) +4*r[1]*r[1]*(Phi[1]+k3[1]*deltaT) -r[2]*r[2]*(Phi[2]+k3[2]*deltaT))/(deltaT*(r[0]*r[0]+0.01*deltaR));
        l4[nR-1] = -0.5*(3*r[nR-1]*r[nR-1]*(Pi[nR-1]+l3[nR-1]*deltaT) -4*r[nR-2]*r[nR-2]*(Pi[nR-2]+l3[nR-2]*deltaT) +r[nR-3]*r[nR-3]*(Pi[nR-3]+l3[nR-3]*deltaT))/(deltaT*r[nR-1]*r[nR-1]);
        for(int ir=1;ir<nR-1;ir++){
            k4[ir] = 0.5*(alpha*(Pi[ir+1]+l3[ir+1]*deltaT)/a -alpha*(Pi[ir-1]+l3[ir-1]*deltaT)/a)/deltaT;
            l4[ir] = 0.5*(r[ir+1]*r[ir+1]*alpha*(Phi[ir+1]+k3[ir+1]*deltaT)/a -r[ir-1]*r[ir-1]*alpha*(Phi[ir-1]+k3[ir-1]*deltaT)/a)/(deltaT*r[ir]*r[ir]);
        }
        //Calculate phi, Phi and Pi on next step
        for(int ir=0;ir<nR;ir++){
            phi[ir] += Pi[ir]*deltaT; //phi is only calculated for visualization purposes
            Phi[ir] += (k1[ir]+2*k2[ir]+2*k3[ir]+k4[ir])*deltaT/6.;
            Pi[ir] += (l1[ir]+2*l2[ir]+2*l3[ir]+l4[ir])*deltaT/6.;
            phi[ir] += Pi[ir]*deltaT+phi[ir]; //phi is only calculated for visualization purposes
        }
    }
    for(int ir=0;ir<(nR/SAVE_RES);ir++){
        Rhistory[ir] = r[ir*SAVE_RES];
    }
    hist[0] = Xhistory;
    hist[1] = Yhistory;
    hist[2] = Rhistory;
    hist[3] = Fhistory;
    return hist;
}


void print_data(double** hist,int print_iterations,int printR){
    FILE* x_data;
    FILE* y_data;
    FILE* r_data;
    FILE* f_data;

    x_data = fopen("Xhistory.dat","w");
    y_data = fopen("Yhistory.dat","w");
    r_data = fopen("Rhistory.dat","w");
    f_data = fopen("Fhistory.dat","w");

    //Print X and Y
    for(int i=0;i<print_iterations-1;i++){
        for(int ir=0;ir<(printR-1);ir++){
            fprintf(x_data,"%lf,",hist[0][i*printR+ir]);
            fprintf(y_data,"%lf,",hist[1][i*printR+ir]);
            fprintf(f_data,"%lf,",hist[3][i*printR+ir]);
        }
        fprintf(x_data,"%lf\n",hist[0][i*printR+printR-1]);
        fprintf(y_data,"%lf\n",hist[1][i*printR+printR-1]);
        fprintf(f_data,"%lf\n",hist[3][i*printR+printR-1]);
    }
    for(int ir=0;ir<(printR-1);ir++){
        fprintf(x_data,"%lf,",hist[0][(print_iterations-1)*printR+ir]);
        fprintf(y_data,"%lf,",hist[1][(print_iterations-1)*printR+ir]);
        fprintf(f_data,"%lf,",hist[3][(print_iterations-1)*printR+ir]);
    }
    fprintf(x_data,"%lf",hist[0][print_iterations*printR-1]);
    fprintf(y_data,"%lf",hist[1][print_iterations*printR-1]);
    fprintf(f_data,"%lf",hist[3][print_iterations*printR-1]);

    //Print R
    for(int ir=0;ir<(printR-1);ir++){
        fprintf(r_data,"%lf,",hist[2][ir]);
    }
    fprintf(r_data,"%lf",hist[2][printR-1]);

    fclose(x_data);
    fclose(y_data);
    fclose(r_data);
    fclose(f_data);

}


void main(){
    double** initial_conditions;
    double** hist;
    double* r;
    double* phi;
    double* Phi;
    double* Pi;
    
    //Define initial conditions
    float p0 = 0.001;
    float r0 = 20.;
    float d = 3.;
    float q = 2.;
    float deltaR = 0.001;
    int maxR = 50;
    initial_conditions = initialize_field(p0,r0,d,q,deltaR,maxR);
    r = initial_conditions[0];
    phi = initial_conditions[1];
    Phi = initial_conditions[2];
    Pi = initial_conditions[3];

    //Pass initial conditions to iteration
    hist = iteration(r,phi,Phi,Pi,deltaR,maxR,ITERATIONS,SAVE_ITERATION);

    //Print simulation history to a file
    int printR = (maxR/deltaR)/SAVE_RES;
    print_data(hist,ITERATIONS/SAVE_ITERATION,printR);

}




