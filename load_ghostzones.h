#pragma once

#include <stdlib.h>
#include <stdio.h>
#include "constants.h"
#include "derivatives.h"


double interp_pol4(double *coarse_fun, int i_coarse){
	int n = 4;
    double x[n], y[n], xp, yp=0, p;

    //Load data from coarser grid
	for(int i=0;i<n;i++){
          x[i] = i+i_coarse-1;
          y[i] = coarse_fun[i+i_coarse-1];
	}

    xp = ((double)i_coarse)+0.5;

	//Calculate Lagrange interpolation
	for(int i=0;i<n;i++){
	    p=1;
		for(int j=0;j<n;j++){
			if(i!=j)
			    p = p*(xp - x[j])/(x[i] - x[j]);
		}
		yp = yp + p*y[i];
	}

    return yp;
}

double interp_pol2(double *coarse_fun, int i_coarse){
    double yp;
    //int i_coarse = (i_fine-1)/GRID_RATIO;
    yp = 0.5*(coarse_fun[i_coarse]+coarse_fun[i_coarse+1]);

    return yp;
}

//Avg overlapping section of the grid with coarse grid (use fine grid_n)
void avg_overlap(double ***all_subgrids,int grid_n){
    double **grid_l = all_subgrids[grid_n];
    double *r = grid_l[0];
    double *phi = grid_l[1];
    double *Phi = grid_l[2];
    double *Pi = grid_l[3];
    double *a = grid_l[4];
    double *alpha = grid_l[5];

    int nR = (int)(*(grid_l[8]));
    bool load_coarser = !((bool)(*(grid_l[12])));
    
    if(load_coarser){
        double     *r_coarse = all_subgrids[grid_n-1][0];
        double   *Phi_coarse = all_subgrids[grid_n-1][2];
        double    *Pi_coarse = all_subgrids[grid_n-1][3];
        double     *a_coarse = all_subgrids[grid_n-1][4];
        double *alpha_coarse = all_subgrids[grid_n-1][5];

        //Ponderate overlapping values from finer grid
        for(int ig=0;ig<OVERLAP;ig++){
            int idx = nR-EXTRA_N + GRID_RATIO*ig;
            //printf("Ponderating r %lf and %lf\n",r[idx],r_coarse[ig]);
            double coarse_weight = (ig+1)/(OVERLAP+1);
            double current_weight = 1-coarse_weight;
            double temp;
            temp = current_weight*Phi[idx] + coarse_weight*Phi_coarse[ig];
            Phi[idx] = temp;
            Phi_coarse[ig] = temp;
            temp = current_weight*Pi[idx] + coarse_weight*Pi_coarse[ig];
            Pi[idx] = temp;
            Pi_coarse[ig] = temp;
            temp = current_weight*a[idx] + coarse_weight*a_coarse[ig];
            a[idx] = temp;
            a_coarse[ig] = temp;
            temp = current_weight*alpha[idx] + coarse_weight*alpha_coarse[ig];
            alpha[idx] = temp;
            alpha_coarse[ig] = temp;
        }
        //Average on the current grid after ponderating
        for(int ig=0;ig<OVERLAP-1;ig++){
            int idx = nR-EXTRA_N + GRID_RATIO*ig;
            //printf("Averaging r %lf and %lf for %lf\n",r[idx],r[idx+GRID_RATIO],r[idx+1]);
              Phi[idx+1] = 0.5*(  Phi[idx] +   Phi[idx+GRID_RATIO]);
               Pi[idx+1] = 0.5*(   Pi[idx] +    Pi[idx+GRID_RATIO]);
                a[idx+1] = 0.5*(    a[idx] +     a[idx+GRID_RATIO]);
            alpha[idx+1] = 0.5*(alpha[idx] + alpha[idx+GRID_RATIO]);
        }
    }
}

void left_ghostzones(double ***all_subgrids,int grid_n){

    double **grid_l = all_subgrids[grid_n];
    double *r = grid_l[0];
    double *phi = grid_l[1];
    double *Phi = grid_l[2];
    double *Pi = grid_l[3];
    double *a = grid_l[4];
    double *alpha = grid_l[5];

    int nR = (int)(*(grid_l[8]));
    bool   load_finer = !((bool)(*(grid_l[11])));

    if(load_finer){
        double     *r_fine = all_subgrids[grid_n+1][0];
        double   *Phi_fine = all_subgrids[grid_n+1][2];
        double    *Pi_fine = all_subgrids[grid_n+1][3];
        double     *a_fine = all_subgrids[grid_n+1][4];
        double *alpha_fine = all_subgrids[grid_n+1][5];
        int        nR_fine = (int)(*(all_subgrids[grid_n+1][8]));

        //Ponderate overlapping values from finer grid
        //for(int ig=0;ig<OVERLAP;ig++){
        //    double fine_weight = (ig+1)/(OVERLAP+1);
        //    double current_weight = 1-fine_weight;
        //    //double current_weight = (ig+1)/(OVERLAP+1);
        //    //double fine_weight = 1-current_weight;
        //      Phi[ig] = current_weight*  Phi[ig] + fine_weight*  Phi_fine[nR_fine-EXTRA_N + GRID_RATIO*ig];
        //       Pi[ig] = current_weight*   Pi[ig] + fine_weight*   Pi_fine[nR_fine-EXTRA_N + GRID_RATIO*ig];
        //        a[ig] = current_weight*    a[ig] + fine_weight*    a_fine[nR_fine-EXTRA_N + GRID_RATIO*ig];
        //    alpha[ig] = current_weight*alpha[ig] + fine_weight*alpha_fine[nR_fine-EXTRA_N + GRID_RATIO*ig];
        //}
        //Get values from finer grid to the ghostzone
        for(int ig=1;ig<=GHOST_SIZE;ig++){
            //printf("copying from r %lf to %lf\n",r_fine[nR_fine-EXTRA_N -GRID_RATIO*ig],r[-ig]);
              Phi[-ig] =   Phi_fine[nR_fine-EXTRA_N -GRID_RATIO*ig];
               Pi[-ig] =    Pi_fine[nR_fine-EXTRA_N -GRID_RATIO*ig];
                a[-ig] =     a_fine[nR_fine-EXTRA_N -GRID_RATIO*ig];
            alpha[-ig] = alpha_fine[nR_fine-EXTRA_N -GRID_RATIO*ig];
        }
    }
}

void right_ghostzones(double ***all_subgrids,int grid_n){

    double **grid_l = all_subgrids[grid_n];
    double *r = grid_l[0];
    double *phi = grid_l[1];
    double *Phi = grid_l[2];
    double *Pi = grid_l[3];
    double *a = grid_l[4];
    double *alpha = grid_l[5];

    int nR = (int)(*(grid_l[8]));
    bool load_coarser = !((bool)(*(grid_l[12])));

    
    if(load_coarser){
        double     *r_coarse = all_subgrids[grid_n-1][0];
        double   *Phi_coarse = all_subgrids[grid_n-1][2];
        double    *Pi_coarse = all_subgrids[grid_n-1][3];
        double     *a_coarse = all_subgrids[grid_n-1][4];
        double *alpha_coarse = all_subgrids[grid_n-1][5];

        //Get values from coarser grid (to the right)
        //  Phi_coarse[0] = 0.5*(  Phi[nR-1] +   Phi_coarse[0]);
        //   Pi_coarse[0] = 0.5*(   Pi[nR-1] +    Pi_coarse[0]);
        //    a_coarse[0] = 0.5*(    a[nR-1] +     a_coarse[0]);
        //alpha_coarse[0] = 0.5*(alpha[nR-1] + alpha_coarse[0]);
        //  Phi_coarse[0] =   Phi[nR-1];
        //   Pi_coarse[0] =    Pi[nR-1];
        //    a_coarse[0] =     a[nR-1];
        //alpha_coarse[0] = alpha[nR-1];
        for(int ig=1;ig<=GHOST_SIZE/GRID_RATIO;ig++){
            //printf("copying from r %lf to %lf\n",r_coarse[ig-1+OVERLAP],r[nR-1 +GRID_RATIO*ig]);
              Phi[nR-1 +GRID_RATIO*ig] =   Phi_coarse[ig-1+OVERLAP];
               Pi[nR-1 +GRID_RATIO*ig] =    Pi_coarse[ig-1+OVERLAP];
                a[nR-1 +GRID_RATIO*ig] =     a_coarse[ig-1+OVERLAP];
            alpha[nR-1 +GRID_RATIO*ig] = alpha_coarse[ig-1+OVERLAP];
        }
        //Get interpolated values from coarser grid
        for(int ig=0;ig<GHOST_SIZE/GRID_RATIO;ig++){
            int idx = nR+GRID_RATIO*ig;
            //printf("copying from r %lf to %lf\n",0.5*(r_coarse[ig-1+OVERLAP]+r_coarse[ig+OVERLAP]),r[idx]);
            //  Phi[idx] = 0.5*(  Phi[idx-1]+  Phi[idx+1]);
            //   Pi[idx] = 0.5*(   Pi[idx-1]+   Pi[idx+1]);
            //    a[idx] = 0.5*(    a[idx-1]+    a[idx+1]);
            //alpha[idx] = 0.5*(alpha[idx-1]+alpha[idx+1]);
              Phi[idx] = interp_pol2(  Phi_coarse, ig-1+OVERLAP);
               Pi[idx] = interp_pol2(   Pi_coarse, ig-1+OVERLAP);
                a[idx] = interp_pol2(    a_coarse, ig-1+OVERLAP);
            alpha[idx] = interp_pol2(alpha_coarse, ig-1+OVERLAP);
        }
    }
}
