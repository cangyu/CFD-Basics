# CFD Basics
Coding practice of Anderson's CFD book: __Computational Fluid Dynamics: The basics with applications__

## Cavity
    TODO

## Laval
The Laval pipe, a classical 1D problem, based on Euler equation.
### 0-Subsonic-Supersonic Isentropic Flow
    MacCormack Scheme is applied.
Usage:
> * Compile: `g++ main.cc -o Laval`
> * Execute: `./Laval`
> * Animate: `python3 animate`

The program will produce a flowfield history file named `flow.txt`, and the steady-stage flowfield looks like:  
![steady-laval](Laval/0/steady.png)
Pay attention to B.C. at both inlet and outlet!

### 1-Subsonic Isentropic Flow
    MacCormack Scheme.  
