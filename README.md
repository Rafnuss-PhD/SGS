# Sequential Gaussian Simualtion

Ever wanted to simulation a Gaussian field easily and fast? This script writtent in MATLAB let you create conditional and unconditional Gaussian realizations.

## Where to start?
1. If you are new to Sequential Gaussian Simulation, it might be a good idea to start by there. I provide a short introduction below with a video illustrative.
2. You should have a look at the Live Script which describe how to use the different function.
3. If you're more interested in what's the difference between the various simulation path and which one is the best adapted for you, have a look at the section below
4. If you wander why is a constant path and what's for.


## What is SGS?
SGS stands for Sequential Gaussian Simulation, as its name suggest, it is a simulation algorithm which generate MultiGaussian field in an iterative manner. Mathematically, it's purpuse is to create a realization ![equation](http://latex.codecogs.com/gif.latex?z%5E%7B%28l%29%7D%28%5Cmathbf%7Bu%7D%29) of a Random Variable ![equation](http://latex.codecogs.com/gif.latex?Z%28%5Cmathbf%7Bu%7D%29%20%5Csim%20%5Cmathcal%7BN%7D%28%5Cboldsymbol%5Cmu_Z%2C%20%5Cboldsymbol%7BC%7D_Z%29). The overall algo can be nicely summerized with:

![equation](http://latex.codecogs.com/gif.latex?Z%20%28%5Cboldsymbol%7Bu%7D_i%29%20%3D%20%5Csum_%7Bj%3D1%7D%5E%7Bi-1%7D%20%5Clambda_j%28%5Cboldsymbol%7Bu%7D_i%29%20Z%28%5Cboldsymbol%7Bu%7D_j%29%20&plus;%20%5Csigma_E%20%28%5Cboldsymbol%7Bu%7D_i%29%20U%28%5Cboldsymbol%7Bu%7D_i%29%2C%20%5Cquad%20%5Cforall%20i%3D1%2C%20%5Cldots%2C%20n%2C)

where ![equation](http://latex.codecogs.com/gif.latex?U) is a standard Gaussian vector responsible to the randomness in the sampling of each value and ![equation](http://latex.codecogs.com/gif.latex?5Clambda_j) are the kriging weightw.


## What's in this package?
This package provides different version/implementation of Sequential Gaussian Simulation, each tailored for another usage.
- [``SGS_trad.m``](https://raphael-nussbaumer-phd.github.io/SGS/html/SGS_trad)
- [``SGS_cst_par.m``](https://raphael-nussbaumer-phd.github.io/SGS/html/SGS_cst_par)
- [``SGS_cst_par_cond.m``](https://raphael-nussbaumer-phd.github.io/SGS/html/SGS_cst_par_cond)
- [``SGS_hybrid.m``](https://raphael-nussbaumer-phd.github.io/SGS/html/SGS_hybrid)



## Which Simulation path to use for SGS?

Nussbaumer, Raphaël, Grégoire Mariethoz, Erwan Gloaguen, and Klaus Holliger. 2017. “Which Path to Choose in Sequential Gaussian Simulation.” _Mathematical Geosciences_. Retrieved (http://link.springer.com/10.1007/s11004-017-9699-5).

## SGS with a Constant Path



## Support or Contact

Having trouble with the code? Contact me.
