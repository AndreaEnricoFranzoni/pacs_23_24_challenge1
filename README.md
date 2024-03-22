#Requirements
GetPot
json
muparser

# Aim
The repo contains a program that wants to find an approximation of the point of minimum for a function that has a 2-dim input and scalar image.
The algos implemented to purse this are: gradient descent with Armijo rule(0), with exponential decay(1), with inverse decay(2), heavy ball(3) and Nesterov(4).
In the first three cases the learning rate is ajusting all along the algorithm, while for the other two is fixed. Sometimes, depending on the parameters' choice, (1) has some difficulties with the second component estimate.

# Parameters
The parameters are passed using a .pot or a .json file:
-f: string that corresponds to the function
-df_x1: string that corresponds to the partial derivative wrt the first variable
-df_x2: string that corresponds to the partial derivative wrt the second variable
NB: since is used muparser to evaluate each one of this three scalar function: write the first variable as 'x1', the second as 'x2'
-x1_0: first component of the initial point (double)
-x2_0: second component of the initial point (double)
-alpha_0: initial value of the learning rate (double) (advice: put it equal to 0.001 for {0,3,4} and equal to 0.1 for {2,3})
-eps_r: tolerance on the residuals (double)
-eps_s: tolerance on the step (double)
-mu: decay parameter for (1) and (2) (double)
-sigma: parameter for (0) (double)
-eta: memory param for (3) and (4) (double)
-itermax: max number of iterations (int)
-rule: which method is used ( 0,1,2,3,4 as above)
-exact_grad: true if we relay on exact gradient, false if we relay on finite differences approximation of it
-h: spacing for finite difference method (double)

# Compilation and execution
'make' compiles everything
'./main' exectues: parameters are by default from .pot file, with -p name_file from another one (.json)
The execution gives back the approximated minimum, along with how quick the convergence has been, or a print message if convergence has not been reached, with a default value of (0.,0.) as result. An error if a wrong algo is passed as parameter (integer different from {0,1,2,3,4}.

# Files
parameters.hpp and ReadParamters.hpp allow to read parameters.
utils.hpp contains some simple but useful generic funcitons (norm evaluation, etc)
challenge1.hpp contains the implementation of the 5 methods, along with a simple way to call them. The methods' names are collected in an enumerate.
in the main.cpp, a big if sequences is done to decide which method is used (a more detailed explanation of why we need it in this way in the file)
