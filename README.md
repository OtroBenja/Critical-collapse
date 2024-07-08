The c files execute the simulation, from a given initial condition

Minkowski.c - simulates on a minkowski spacetime with Euler integration

MinkowskiRK4.c - simulates on a minkowski spacetime with RK4 integration

RK4.c - simulates on the spacetime given by Choptuik (1993) with RK4 integration



Every simulation outputs 3 files:

Xhistory.dat - A history of the scalar field X, a transformation of Φ

Yhistory.dat - A history of the scalar field Y, a transformation of Π

Rhistory.dat - The values of the radial space coordinate, r



Not every time step is saved, nor every value of r. 
This is to save time during the execution and the rendering of the animation,
as well as reducing the file size of the output. How any values are saved can be
changed with the constants SAVE_RES and SAVE_ITERATION in the C file.


The file plot.py loads the previous output and renders an animation of the simulation.

