*********LSN_Exercises_ Delivery*********
+ @author Giovanni Ferrami
+ @e-mail g.ferrami@gmail.com

Each exercise has a dedicated folder, inside it open "LSN_Solution_**.ipynb" to see the commented solution of the assigned exercises (the exercise instructions can be found in "LSN_Exercise_**.ipynb").
Each exercise is implemented in C++ and then the generated data is analyzed in Python3 through the notebooks uploaded here.


Here's a list of the exercises with a brief comment on what they are about:
- Ex01: 
        Generating simple statistical distribution true the inverse of the cumulative function. 
        Testing the Central Limit Theorem and stable distributions.
        Simulation of Buffon needle experiment to give a rough estimate of Pi.
- Ex02: 
        Uniform sampling vs Importance sampling
        Discrete and continous 3D-random walks mean distance from the origin.
- Ex03: 
        Simple finance simulation. Estimate the Put and Call value of an European option, comparing it with the analitic result from Black and Scholes algorithm
- Ex04: 
        Molecular Dynamics simulation with Verlet algorithm of a set of 100+ molecules interacting through a Lennard-Jones potential in periodic      boundary conditions.
        Implementing a thermal stabilization algorithm, step by step evaluation of kinetic and potential energies, temperature, pressure.
        Applying the simulation to Argon and Krypton in solid, liquid and gaseous phases, retrieving coherent phisical values for the observed quantities.
- Ex05: 
        Metropolis algorithm. Evaluating hydrogenic probability distribution for 1S and 2P energy levels.
- Ex06: 
        1-D Ising model simulated through Metropolis and Gibbs sampling algorithms. Evaluating some termodynamic parameter confronting the two algos with the analitic solution.
- Ex07: 
        Monte Carlo NVT applied to the same Molecular Dynamics system analyzed in ex04.
- Ex08: 
        Variational Monte Carlo code for a single quantum particle in 1D which exploits the Metropolis algorithm to sample a trial wave function $|\Psi_T(x)|^2$
- Ex09: 
        Solving the Traveling Salesman Problem (TSP) through a genetic algorithm written from scratch in C++
- Ex10: 
        Adapting the Genetic Algorithm code, developed during the Numerical Exercise 9, to solve the TSP with a Simulated Annealing (SA) algorithm.
        In the second part we converted the SA algorithm from serial to parallel computing, testing the performances of a multi-CPU code.
- Ex11:
        Performing machine learning regression on noisy data with a Neural Network (NN) (*only python)
- Ex12:
        Using deep convolutional neural network models, implemented in the Keras python package, to recognize and distinguish between the ten handwritten digits (0-9) (*only python)