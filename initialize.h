#include "derivatives.h"
#include "constants.h"


double** initialize_field(int fType,double* model_parameters,double deltaR,double minR,double maxR){
    double p0 = model_parameters[0];
    double r0 = model_parameters[1];
    double d  = model_parameters[2];
    int nR = (int)((maxR-minR)/deltaR);
    double* r = malloc(sizeof(double)*nR);
    double* phi = malloc(sizeof(double)*nR);
    double* Phi = malloc(sizeof(double)*nR);
    double* Pi = malloc(sizeof(double)*nR);
    double** return_values = malloc(sizeof(double*)*4);

    //Define r and calculate initial phi
    for(int i=0;i<nR;i++)
        r[i] = minR + i*deltaR;
    if(minR==0.0) r[0] = 1.0E-50;

    if(fType==0){
        for(int i=0;i<nR;i++)
            phi[i] = p0*tanh((r[i]-r0)/d);
    }
    if(fType==1){
        for(int i=0;i<nR;i++)
            phi[i] = p0*pow(r[i],3)*pow(E,-pow((r[i]-r0)/d,2));
    }
    if(fType==2){
        for(int i=0;i<nR;i++)
            phi[i] = p0;
    }
    if(fType==3){
        for(int i=0;i<nR;i++)
            phi[i] = p0*r[i];
    }

    //Calculate initial Phi
    Phi[0] = leftmost_D1(phi, 0, deltaR);
    Phi[1] = leftmid_D1(phi, 1, deltaR);
    for(int i=2;i<nR-2;i++){
        Phi[i] = centered_D1(phi, i, deltaR);
    }
    Phi[nR-2] = rightmid_D1(phi, nR-2, deltaR);
    Phi[nR-1] = rightmost_D1(phi, nR-2, deltaR);

    //Set initial Pi to zero
    for(int i=0;i<nR;i++){
        Pi[i] = 0;
    }

    
    //Set initial Pi to move left
    //double* rphi = malloc(sizeof(double)*nR);
    //for(int i=0;i<nR;i++){
    //    rphi[i] = r[i]*phi[i];
    //}
    //if(minR==0.0) {Pi[0] = 0;}
    //else Pi[0] = leftmost_D1(rphi, 0, deltaR)/r[0];
    //Pi[1] = leftmid_D1(rphi, 1, deltaR)/r[1];
    //for(int i=2;i<nR-2;i++){
    //    Pi[i] = centered_D1(rphi, i, deltaR)/r[i];
    //}
    //Pi[nR-2] = rightmid_D1(rphi, nR-2, deltaR)/r[nR-2];
    //Pi[nR-1] = rightmost_D1(rphi, nR-1, deltaR)/r[nR-1];
    //free(rphi);
        
    return_values[0] = r;
    return_values[1] = phi;
    return_values[2] = Phi;
    return_values[3] = Pi;
 
    return return_values;
}