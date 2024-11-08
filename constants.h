#pragma once

#define MASS 1
#define METRIC 1 // 0 = minkowski; 1 = choptuik; 2 = modified choptuik
#define SAVE_MODE 0 // 0 = save uniformly on every SAVE_RES and SAVE_ITERATION ; 1 = save all points after FIRST_ITERATION and with r > MIN_R
#define SAVE_RES 50
#define SAVE_ITERATION 100
#define FIRST_ITERATION 77000 //77950
#define MIN_R 0
#define ITERATIONS 78500
#define EPSILON 0.0 // Default Kreiss-Oliger dampening
#define BETA 0.0
#define PI 3.141592653
#define _6 0.166666666 // 1/6
#define E  2.718281828