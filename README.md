The C files execute the simulation, from a given initial condition

- Colapso_RK4.c - uses RK4 integration to simulate on the spacetime given by Choptuik (1993) or minkowski, this behaviour can
be modified by setting the constant MINKOWSKI to 0 (false) or 1 (true) before compilation.

The program can receive some arguments when being run and they are all optional, if the value given is 0 then the default 
value is assigned. The arguments are the following:

- fType      (int)  Type of the initial function for $\phi$, can be either 0 (tanh) or 1 (gaussian)
- p0       (float)  Amplitude of the initial function
- r0       (float)
- d        (float)
- q        (float)
- maxR       (int)  Maximun value for the radial coordinate, determinates the total space to simulate on
- iterations (int)  Total number of iterations to run the simulation
- nThreads   (int)  Number of threads used for paralelization

The program outputs 2 files after finishing, the first one contains the following:
- Various information about the parameters given and the runtime
- The values of the radial space coordinate, r
- A history of the scalar field φ
- A history of the scalar field Φ, the derivative dφ/dr
- A history of the scalar field Π, the derivative (a/alpha)*(dφ/dt)

The second one contains the following:
- The mass of the field in every saved iteration, this is supposed to be a constant


Not every time step is saved, nor every value of r. 
This is to save time during the execution and the rendering of the animation,
as well as reducing the file size of the output. How many values are saved can be
changed with the constants SAVE_RES and SAVE_ITERATION in the C file.


The notebook plot.ipynb can be used to load previous outputs and see a speed comparison as well as render an animation of the simulation

