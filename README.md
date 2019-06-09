# CFD Basics
Coding practice of Anderson's CFD book: __Computational Fluid Dynamics: The basics with applications__

## Laval
    The Laval pipe, a classical 1D problem, based on Euler equation.
### 0-Subsonic-Supersonic Isentropic Flow
    MacCormack Scheme.
Usage:
> * Compile: `g++ main.cc -o Laval`
> * Execute: `./Laval`
> * Animate: `python3 animate`

The program will produce a flowfield history file named `flow.txt`, and the steady-state flowfield looks like:  
![steady-laval](Laval/0/steady.png)
Pay attention to B.C. at both inlet and outlet!

### 1-Subsonic Isentropic Flow
    MacCormack Scheme.  
    Clearly, velocity peaks at central.
    
Usage:
> * Compile: `g++ main.cc -o Laval`
> * Execute: `./Laval`
> * Animate: `python3 animate`

### 2-Conservative form for Subsonic-Supersonic Isentropic Flow
    MacCormack Scheme.
    
Usage:
> * Compile: `g++ main.cc -o Laval`
> * Execute: `./Laval`
> * Animate: `python3 animate`

### 3-Shockwave Capture
    MacCormack Scheme.  
	Add artificial viscosity at both prediction and correction steps.

Usage:
> * Compile: `g++ main.cc -o Laval`
> * Execute: `./Laval`
> * Animate: `python3 animate`

The program will produce a flowfield history file named `flow.txt`, and the steady-state flowfield looks like:  
![steady-shock](Laval/3/steady.png)

## Cavity
    TODO