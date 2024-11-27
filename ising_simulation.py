#Translated julia script to python
#In julia you can use greek symbols as variables and some other operators, you'll probably have to parse all these symbols to strings
import numpy as np
from numba import njit #Maybe you don't need this package, i was using multithreading to simulate one sequence on each core, there is probably a better way to do this in python

#There is probably a better way to do an ising simulation with Glauber dynamics in python

# Simulates the Ising model in 2D using Glauber dynamics.

# Generates a random 2D Ising configuration with dimensions m x n.
# Elements are randomly assigned -1 or 1.
def rand_ising2d(m, n=None):
    if n is None:
        n = m
    return np.random.choice([-1, 1], size=(m, n), dtype=np.int8)

# Set random seed for reproducibility.
np.random.seed(4649)

# Critical inverse temperature for the 2D Ising model.
# β (beta) is replaced with beta_crit
beta_crit = np.log(1 + np.sqrt(2)) / 2

# Precompute probabilities for spin flips at the critical temperature.
# This avoids repeated calculations within the evolution function.
prob = np.exp(-2 * beta_crit * np.arange(-4, 5))

# Evolves a 2D Ising configuration `s` for `niters` iterations at inverse temperature `beta`.
# β (beta) is replaced with beta
@njit  # Use Numba for JIT compilation and speedup
def ising2d_evolve(s, beta, niters):
    m, n = s.shape  # Get dimensions of the lattice
    prob = np.exp(-2 * beta * np.arange(-4, 5))  # Probabilities for spin flips
    for _ in range(niters):  # Loop over iterations
        for j in range(n):
            for i in range(m):
                # Get the neighboring spins with periodic boundary conditions
                NN = s[(i - 1) % m, j]
                SS = s[(i + 1) % m, j]
                WW = s[i, (j - 1) % n]
                EE = s[i, (j + 1) % n]
                CT = s[i, j]
                # Calculate the energy change if the spin is flipped
                k = CT * (NN + SS + WW + EE)
                # Flip the spin with probability prob[k+4] (adjusting index for Python)
                if np.random.rand() < prob[k + 4]:
                    s[i, j] *= -1
    return s

# Generates an initial ensemble of M random Ising configurations with dimensions N x N.
def isingInitialEnsemble(M, N):
    return [rand_ising2d(N) for _ in range(M)]

# Generates a dataset of M Ising configurations, each evolving for n_frames frames.
# Starts with an initial random configuration and quenches to inverse temperature `beta`.
# β (beta) is replaced with beta
def isingDataSet(M, N, beta, n_frames, n_t=1):
    X = []  # Initialize an empty list to store configurations
    X0 = isingInitialEnsemble(M, N)  # Generate initial random configurations

    for i in range(M):
        X.append(X0[i].copy())  # Add the initial configuration to the dataset
        for _ in range(1, n_frames):
            X.append(X[-1].copy())  # Add a copy of the previous configuration
            X[-1] = ising2d_evolve(X[-1], beta, n_t)  # Evolve the configuration for n_t steps

    return X

# Computes the distance matrix between Ising configurations in A.
# The distance between two configurations is the sum of squared differences between their spins.
def isingDistanceMatrix(A):
    M = len(A)
    D = np.zeros((M, M))
    for i in range(M):
        for j in range(i):  # Python uses 0-based indexing
            D[i, j] = np.sum((A[i] - A[j]) ** 2)  # Calculate the distance
            D[j, i] = D[i, j]  # Ensure symmetry
    return D

# Computes the Diffusion Map distance matrix between Ising configurations in A.
# Uses the exponential of the negative squared distance scaled by epsilon.
# ϵ (epsilon) is replaced with epsilon
def isingDiffMapDistance(A, M, N, epsilon):
    D = np.zeros((M, M))  # Initialize the distance matrix
    for i in range(M):
        for j in range(M):
            # Calculate the Diffusion Map distance
            D[i, j] = np.exp(-isingDistance(A[i], A[j]) / ((2 * N) ** 2 * epsilon))  
    return D

# Generates a dataset of M Ising configurations, each quenched to inverse temperature `beta`.
# Only returns the final state of each realization.
# β (beta) is replaced with beta
def isingSteadyDataSet(M, N, beta, n_t):
    X = isingInitialEnsemble(M, N)  # Generate initial random configurations
    for x in X:
        x = ising2d_evolve(x, beta, n_t)  # Evolve the configuration for n_t steps
    return X