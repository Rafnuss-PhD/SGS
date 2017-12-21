# Sequential Gaussian Simulation

Have you ever wanted to generate a Gaussian field? This MATLAB script let you to easily create multiple conditional or unconditional 2D realizations of a Gaussian model of your choice.



## Where to start?
The LiveScript [``script_SGS.m``](https://raphael-nussbaumer-phd.github.io/SGS/html/script_SGS) gives a nice overview of the different implementation of SGS available in this code. 

Code is structure as follows. [``SGS.m``](https://raphael-nussbaumer-phd.github.io/SGS/html/SGS) is the main function, it explains the input and outpur format, it also check the input argument and redirect toward the appropriate implementation of SGS, as listed below:

This package provides different implementation of Sequential Gaussian Simulation, each tailored for another usage.
- [``SGS_trad.m``](https://raphael-nussbaumer-phd.github.io/SGS/html/SGS_trad): the traditional alogrithm.
- [``SGS_cst.m``](https://raphael-nussbaumer-phd.github.io/SGS/html/SGS_cst): the constant path algorithm.
- [``SGS_cst_par.m``](https://raphael-nussbaumer-phd.github.io/SGS/html/SGS_cst_par): the constant path algorithm with parralelisation. 
- [``SGS_cst_par_cond.m``](https://raphael-nussbaumer-phd.github.io/SGS/html/SGS_cst_par_cond): the constant path algorithm with conditional point.
- [``SGS_hybrid.m``](https://raphael-nussbaumer-phd.github.io/SGS/html/SGS_hybrid): Hydrid constant/randomized path switching at a certain level of the Multi-grid path.
- [``SGS_varcovar.m``](https://raphael-nussbaumer-phd.github.io/SGS/html/SGS_varcovar): Compute the full covariance matrice of the simulation.



## What is SGS?
SGS stands for Sequential Gaussian Simulation, as its name suggest, it is a simulation algorithm which generate MultiGaussian field in an iterative manner. Mathematically written, it's purpuse is to create a realization ![equation](http://latex.codecogs.com/gif.latex?z%5E%7B%28l%29%7D%28%5Cmathbf%7Bu%7D%29) from a random variable ![equation](http://latex.codecogs.com/gif.latex?Z%28%5Cmathbf%7Bu%7D%29%20%5Csim%20%5Cmathcal%7BN%7D%28%5Cboldsymbol%5Cmu_Z%2C%20%5Cboldsymbol%7BC%7D_Z%29). The overall algo can be nicely summerized with:

![equation](http://latex.codecogs.com/gif.latex?Z%20%28%5Cboldsymbol%7Bu%7D_i%29%20%3D%20%5Csum_%7Bj%3D1%7D%5E%7Bi-1%7D%20%5Clambda_j%28%5Cboldsymbol%7Bu%7D_i%29%20Z%28%5Cboldsymbol%7Bu%7D_j%29%20&plus;%20%5Csigma_E%20%28%5Cboldsymbol%7Bu%7D_i%29%20U%28%5Cboldsymbol%7Bu%7D_i%29%2C%20%5Cquad%20%5Cforall%20i%3D1%2C%20%5Cldots%2C%20n%2C)

where ![equation](http://latex.codecogs.com/gif.latex?U) is a standard Gaussian vector responsible to the randomness in the sampling of each value and ![equation](http://latex.codecogs.com/gif.latex?%5Clambda_j) are the kriging weights.



## Which Simulation path to use for SGS?
My first study of SGS focused on better understanding the effect of the different simulation path, that is, the order in which the nodes are simulated. We also tried to provide some guideline on which path to avoir and which one to favor. In a nutshell, a simulation which maximizes the overall distance among the nodes during the simulation shows better performance than the one simulating consecutive nodes. 

Have a look at the paper if you're interested, otherwise, have a look at [a presentation I gave at Geostatistics Valencia 2016](https://www.researchgate.net/publication/318858970_Sequential_Simulation_Path_Biases_and_how_to_live_with_them)

Nussbaumer, Raphaël, Grégoire Mariethoz, Erwan Gloaguen, and Klaus Holliger. 2017. “Which Path to Choose in Sequential Gaussian Simulation.” _Mathematical Geosciences_. Retrieved DOI:[10.1007/s11004-017-9699-5](http://link.springer.com/10.1007/s11004-017-9699-5).


## SGS with a Constant Path
In this second paper, I looked at a quite oftenly used technique which consist of using the same simulation path among multiple realizations. This is done in order to be able to reuse the same kriging weights for all realization, which in turn allows a quite drastic reduction of the computational cost. In this work we wanted to quantify the amount of bias added with this technique, and provide some recommendation on the implementation. We showed that these biases were quite unsignificant compared to the computational gain obtained. Have a look a the paper below for more detail.
The script [``cst_path_paper/script_paper.m``](https://raphael-nussbaumer-phd.github.io/SGS/cst_path_paper/html/script_paper) shows you how the figure were built.

Nussbaumer, Raphaël, Grégoire Mariethoz, Mathieu Gravey, Erwan Gloaguen, and Klaus Holliger. 2017. “Accelerating Sequential Gaussian Simulation with a Constant Path.” Computers & Geosciences. Retrieved DOI:[10.1016/j.cageo.2017.12.006](http://linkinghub.elsevier.com/retrieve/pii/S0098300417304685).


## Support or Contact
Having trouble with the code? Contact me.

More publication and presentation on researchgate:
[![Image](https://upload.wikimedia.org/wikipedia/commons/a/aa/ResearchGate_Logo.png)](https://www.researchgate.net/project/Simulation-Path-in-Sequential-Gaussian-Simulation)
