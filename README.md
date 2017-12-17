# Sequential Gaussian Simulation

Have you ever wanted to generate a Gaussian field? This MATLAB script let you to easily create multiple conditional or unconditional 2D realizations.


## Where to start?
1. If you are new to Sequential Gaussian Simulation, it might be a good idea to start by there. I provide a short introduction below with a video illustrative.
2. You should have a look at the Live Script which describe how to use the different function.
3. If you're more interested in what's the difference between the various simulation path and which one is the best adapted for you, have a look at the section below
4. If you wander why is a constant path and what's for.


## What is SGS?
SGS stands for Sequential Gaussian Simulation, as its name suggest, it is a simulation algorithm which generate MultiGaussian field in an iterative manner. Mathematically written, it's purpuse is to create a realization ![equation](http://latex.codecogs.com/gif.latex?z%5E%7B%28l%29%7D%28%5Cmathbf%7Bu%7D%29) from a random variable ![equation](http://latex.codecogs.com/gif.latex?Z%28%5Cmathbf%7Bu%7D%29%20%5Csim%20%5Cmathcal%7BN%7D%28%5Cboldsymbol%5Cmu_Z%2C%20%5Cboldsymbol%7BC%7D_Z%29). The overall algo can be nicely summerized with:

![equation](http://latex.codecogs.com/gif.latex?Z%20%28%5Cboldsymbol%7Bu%7D_i%29%20%3D%20%5Csum_%7Bj%3D1%7D%5E%7Bi-1%7D%20%5Clambda_j%28%5Cboldsymbol%7Bu%7D_i%29%20Z%28%5Cboldsymbol%7Bu%7D_j%29%20&plus;%20%5Csigma_E%20%28%5Cboldsymbol%7Bu%7D_i%29%20U%28%5Cboldsymbol%7Bu%7D_i%29%2C%20%5Cquad%20%5Cforall%20i%3D1%2C%20%5Cldots%2C%20n%2C)

where ![equation](http://latex.codecogs.com/gif.latex?U) is a standard Gaussian vector responsible to the randomness in the sampling of each value and ![equation](http://latex.codecogs.com/gif.latex?%5Clambda_j) are the kriging weights.


## What's in this package?
This package provides different implementation of Sequential Gaussian Simulation, each tailored for another usage.
- [``SGS_trad.m``](https://raphael-nussbaumer-phd.github.io/SGS/html/SGS_trad): the traditional alogrithm.
- [``SGS_cst_par.m``](https://raphael-nussbaumer-phd.github.io/SGS/html/SGS_cst_par): the constant path algorithm with parralelisation. 
- [``SGS_cst_par_cond.m``](https://raphael-nussbaumer-phd.github.io/SGS/html/SGS_cst_par_cond): the constant path algorithm with conditional point.
- [``SGS_hybrid.m``](https://raphael-nussbaumer-phd.github.io/SGS/html/SGS_hybrid): Hydrid constant/randomized path switching at a certain level of the Multi-grid path.


## Which Simulation path to use for SGS?
My first study of SGS focused on better understanding the effect of the different simulation path, that is, the order in which the nodes are simulated. We also tried to provide some guideline on which path to avoir and which one to favor. In a nutshell, a simulation which maximizes the overall distance among the nodes during the simulation shows better performance than the one simulating consecutive nodes. 

Have a look at the paper if you're interested, otherwise, have a look at [a presentation I gave at Geostatistics Valencia 2016](https://www.researchgate.net/publication/318858970_Sequential_Simulation_Path_Biases_and_how_to_live_with_them)

Nussbaumer, Raphaël, Grégoire Mariethoz, Erwan Gloaguen, and Klaus Holliger. 2017. “Which Path to Choose in Sequential Gaussian Simulation.” _Mathematical Geosciences_. Retrieved (http://link.springer.com/10.1007/s11004-017-9699-5).


## SGS with a Constant Path
In this second paper, I looked at the quite oftenly used technique which consist of using the same simulationp path among multiple realizations in order to be able to reuse the same kriging weights for all of them. In this work we wanted to quantify the amount of bias added with this technique, and provide some recommendation on the implementation. We showed that these biases were quite unsignificant compared to the computational gain obtained. 

The paper will soon be available...


## Support or Contact
Having trouble with the code? Contact me.

More publication and presentation on researchgate:
[![Image](https://upload.wikimedia.org/wikipedia/commons/a/aa/ResearchGate_Logo.png)](https://www.researchgate.net/project/Simulation-Path-in-Sequential-Gaussian-Simulation)
