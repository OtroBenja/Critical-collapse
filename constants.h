#pragma once

#define MASS 0 // 1 = Make a mass history file; 0 = Don't
#define METRIC 1 // 0 = minkowski; 1 = choptuik
#define SAVE_MODE 0 // 0 = save uniformly on every SAVE_RES and SAVE_ITERATION ; 1 = save all points after FIRST_ITERATION and with r > MIN_R
#define SAVE_RES 20
#define SAVE_ITERATION 500
#define FIRST_ITERATION 0 //77950
#define MIN_R 30
#define ITERATIONS 7850
#define EPSILON 0.0 // Default Kreiss-Oliger dampening
#define BETA 0.0
#define PI 3.141592653
#define _6 0.166666666 // 1/6
#define E  2.718281828
#define TOLERANCE 1.0E-2 // Min value of g^rr to consider as non-zero for the purposes of a collapse into BH
#define N_LEVELS 2 // Maximum amount of grid levels in variable grid
#define GRID_RATIO 2 // Ratio of resolution from one grid to the next