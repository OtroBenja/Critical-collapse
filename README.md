The C files execute the simulation, from a given initial condition

- MinkowskiRK4.c - simulates on a minkowski spacetime with RK4 integration
- RK4.c - simulates on the spacetime given by Choptuik (1993) with RK4 integration

The program can recive some arguments being run and they are all optional, if the value given is 0 then the default 
value is assigned. The arguments are the following:

- fType      (int)  Type of the initial function for $\phi$, can be either 0 (tanh) or 1 (gaussian)
- p0       (float)  Amplitude of the initial function
- r0       (float)
- d        (float)
- q        (float)
- iterations (int)  Total number of iterations to run the simulation
- maxR       (int)  Maximun value for the radial coordinate, determinates the total space to simulate on
- nThreads   (int)  Number of threads used for paralelization


Every simulation outputs 4 files (This is going to be changed soon to 1 file only):

- Fhistory.dat - A history of the scalar field φ
- Xhistory.dat - A history of the scalar field Φ, the derivative dφ/dr
- Yhistory.dat - A history of the scalar field Π, the derivative (a/alpha)*(dφ/dt)
- Rhistory.dat - The values of the radial space coordinate, r



Not every time step is saved, nor every value of r. 
This is to save time during the execution and the rendering of the animation,
as well as reducing the file size of the output. How many values are saved can be
changed with the constants SAVE_RES and SAVE_ITERATION in the C file.


The file plot.py loads the previous output and renders an animation of the simulation.

