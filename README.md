The C files execute the simulation, from a given initial condition

Minkowski.c - simulates on a minkowski spacetime with Euler integration

MinkowskiRK4.c - simulates on a minkowski spacetime with RK4 integration

RK4.c - simulates on the spacetime given by Choptuik (1993) with RK4 integration

Both RK4 methods have and outgoing wave condition on the outer frontier r=50


Every simulation outputs 4 files:

Fhistory.dat - A history of the scalar field φ

Xhistory.dat - A history of the scalar field Φ, the derivative dφ/dr

Yhistory.dat - A history of the scalar field Π, the derivative (a/alpha)*(dφ/dt)

Rhistory.dat - The values of the radial space coordinate, r



Not every time step is saved, nor every value of r. 
This is to save time during the execution and the rendering of the animation,
as well as reducing the file size of the output. How many values are saved can be
changed with the constants SAVE_RES and SAVE_ITERATION in the C file.


The file plot.py loads the previous output and renders an animation of the simulation.

