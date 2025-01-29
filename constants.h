#pragma once

#define MASS 0 // 1 = Make a mass history file; 0 = Don't
#define METRIC 0 // 0 = minkowski; 1 = choptuik; 2 = quasi-static choptuik
#define SAVE_MODE 0 // 0 = save uniformly on every SAVE_RES and SAVE_ITERATION ; 1 = save all points after FIRST_ITERATION and with r > MIN_R
#define SAVE_RES 20
#define SAVE_ITERATION 10
#define FIRST_ITERATION 8000 //77950
#define MIN_R 48
#define MAX_R 50
#define ITERATIONS 10000 //Max iterations by default
#define BETA 0.0
#define PI 3.141592653 // π number
#define _6 0.166666666 // 1/6
#define E  2.718281828 // e number
#define TOLERANCE 1.0E-2 // Min value of g^rr to consider as non-zero for the purposes of a collapse into BH
#define SUBGRID_MODE 0 // 0 = Subgrid size decreaces exponentially towards the origin ; 1 = All subgrids are the same size
#define N_LEVELS 4 // Maximum amount of grid levels in variable grid
#define GRID_RATIO 2 // Ratio of resolution from one grid to the next (Only works for 2)
#define GHOST_SIZE 8*GRID_RATIO // Size of the ghostzone for each subgrid
#define OVERLAP 1 // Amount of overlapping points between subgrids (min 1)
#define EXTRA_N (1+GRID_RATIO*(OVERLAP-1)) // Amount of extra points for overlapping subgrids (DO NOT CHANGE)