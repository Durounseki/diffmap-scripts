#Vicsek simulation using intrinsic noise in reorientation updating rule and nearest neighbor simulation using a cell list.
#Probably there is an easier and more direct way to do this in python

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# Calculates the nearest neighbors cells in a d-dimensional space.
# Returns a list of tuples representing the neighboring cells.
def neighbors(d: int):
    indices = np.ndindex((3,) * d)  # Create indices for all neighbors including self
    return list(indices)[: len(indices) // 2]  # Return only the first half to avoid double counting

# Partitions the space into cells of size r_0 in a d-dimensional space with side length L.
# Returns a dictionary where keys are tuples representing cell positions and
# values are empty lists to store particle indices.
def CellGrid(d: int, L, r_0):
    indices = np.ndindex((int(L // r_0) + 1,) * d)  # Create indices for all cells including ghost cells
    grid = {index: [] for index in indices}  # Initialize a dictionary to store the grid
    return grid

# Creates a cell list for particles P in a d-dimensional space with side length L and interaction radius r_0.
# Assigns each particle to its corresponding cell in the grid `cells`.
# Includes ghost cells to handle periodic boundary conditions.
def CellList(P, L, r_0, cells):
    n = len(P)  # Number of particles
    for i in range(n):
        # Assign particle i to its cell based on its position
        cell_index = tuple(int(P[i][j] // r_0) for j in range(len(P[i])))
        cells[cell_index].append(i)
    periodicBC(cells, L, r_0)  # Apply periodic boundary conditions to the cell list

# Returns the tuples representing the ghost cells for a d-dimensional space with side length L and cell size r_0.
# Ghost cells are copies of boundary cells used to handle periodic boundary conditions.
def ghostCells(L, r_0):
    max_index = int(L // r_0)
    tbBoundary = [(i, max_index) for i in range(max_index + 1)]  # Top boundary cells
    tbBoundary += [(i, 0) for i in range(max_index + 1)]  # Bottom boundary cells
    rlBoundary = [(max_index, i) for i in range(max_index + 1)]  # Right boundary cells
    rlBoundary += [(0, i) for i in range(max_index + 1)]  # Left boundary cells
    return tbBoundary, rlBoundary

# Applies periodic boundary conditions to the cell list `cells`.
# Copies the particles in boundary cells to their corresponding ghost cells.
def periodicBC(cells, L, r_0):
    max_index = int(L // r_0)
    tbBoundary, rlBoundary = ghostCells(L, r_0)  # Get the ghost cell indices

    # Handle top and bottom boundary conditions
    for boundaryCell in tbBoundary:
        if boundaryCell in cells:  # Check if the boundary cell exists
            if boundaryCell[1] == 0:  # If it's a bottom boundary cell
                cells[(boundaryCell[0], max_index)] = cells[boundaryCell]  # Copy to top ghost cell
            elif boundaryCell[1] == max_index:  # If it's a top boundary cell
                cells[(boundaryCell[0], 0)] = cells[boundaryCell]  # Copy to bottom ghost cell

    # Handle right and left boundary conditions
    for boundaryCell in rlBoundary:
        if boundaryCell in cells:  # Check if the boundary cell exists
            if boundaryCell[0] == 0:  # If it's a left boundary cell
                cells[(max_index, boundaryCell[1])] = cells[boundaryCell]  # Copy to right ghost cell
            elif boundaryCell[0] == max_index:  # If it's a right boundary cell
                cells[(0, boundaryCell[1])] = cells[boundaryCell]  # Copy to left ghost cell

# Calculates the nearest neighbors within a cell or between two cells.
# Updates the `ps` matrix (fixed radius neighbors matrix) based on distances between particles.
def cellNeighbors(ps, Is, p, r_0):
    for k, i in enumerate(Is[:-1]):  # Iterate over all pairs of particles in the cell
        for j in Is[k + 1:]:
            if np.linalg.norm(np.array(p[i]) - np.array(p[j])) <= r_0:  # Check if the distance is within the interaction radius
                ps[i][j] = 1  # Mark as neighbors
                ps[j][i] = 1
            else:
                ps[i][j] = 0
                ps[j][i] = 0
        ps[i][i] = 1  # Each particle is a neighbor of itself

# Calculates the nearest neighbors between two cells.
# Updates the `ps` matrix (fixed radius neighbors matrix) based on distances between particles.
def cellNeighbors2(ps, Is, js, p, r_0):
    for i in Is:  # Iterate over particles in the first cell
        for j in js:  # Iterate over particles in the second cell
            if np.linalg.norm(np.array(p[i]) - np.array(p[j])) <= r_0:  # Check if the distance is within the interaction radius
                ps[i][j] = 1  # Mark as neighbors
                ps[j][i] = 1
            else:
                ps[i][j] = 0
                ps[j][i] = 0
        ps[i][i] = 1  # Each particle is a neighbor of itself

# Computes the nearest neighbors matrix `FRN` for particles P using a cell list.
# FRN_{ij} = 1 if the distance between p_i and p_j is less than the interaction radius r_0.
def near_neighbors(cells, P, r_0, FRN):
    offsets = neighbors(len(P[0]))  # Get the neighboring cell offsets

    # Iterate over non-empty cells
    for cell, Is in cells.items():
        # Pairs of points within the cell
        cellNeighbors(FRN, Is, P, r_0)  # Calculate neighbors within the cell

        # Pairs of points with non-empty neighboring cells
        for offset in offsets:
            neigh_cell = tuple(np.array(cell) + np.array(offset))  # Get the neighboring cell
            if neigh_cell in cells:  # Check if the neighboring cell exists
                js = cells[neigh_cell]  # Get the particles in the neighboring cell
                cellNeighbors2(FRN, Is, js, P, r_0)  # Calculate neighbors between the cells

# Empties all cells in the cell list `cells`.
def emptyCells(cells):
    for cell in cells.values():
        cell.clear()  # Clear the list of particles in each cell

# Evolves the Vicsek model for one time step.
# Updates the positions (X), orientations (Theta), and velocities (V) of the particles.
def fly(X, Theta, V, dt, L, r_0, N, v_0, eta, cells, fixedRN):
    CellList(X, L, r_0, cells)  # Create the cell list
    near_neighbors(cells, X, r_0, fixedRN)  # Compute the nearest neighbors matrix

    # Calculate the mean velocity of neighbors for each particle
    Vmean = [
        np.mean(np.array([V[j] for j in np.where(fixedRN[i] > 0)[0]]), axis=0)
        for i in range(N)
    ]
    Vmean = [v / np.linalg.norm(v) if np.linalg.norm(v) > 0 else np.array([1, 0]) for v in Vmean]  # Normalize

    # Update the orientations with intrinsic noise
    Theta_new = [np.arctan2(u[1], u[0]) + eta * (np.random.rand() - 0.5) for u in Vmean]

    # Update the velocities based on the new orientations
    V_new = np.array([[v_0 * np.cos(tt), v_0 * np.sin(tt)] for tt in Theta_new])

    # Update the positions, orientations, and velocities of the particles
    for i in range(N):
        V[i] = V_new[i]
        Theta[i] = Theta_new[i]
        X[i] = X[i] + dt * V_new[i]

    PosPeriodicBC(X, L)  # Apply periodic boundary conditions to the positions

# Applies periodic boundary conditions to the positions of the particles.
def PosPeriodicBC(X, L):
    for i in range(len(X)):
        X[i][0] = X[i][0] % L  # Apply periodic boundary condition in the x-dimension
        X[i][1] = X[i][1] % L  # Apply periodic boundary condition in the y-dimension

# Simulate one trajectory
def birdTrajectory(X, V, dt, N, nT, r_0, v_0, eta, L, d):
    cells = CellGrid(d, L, r_0)  # Create the cell grid
    fixedRN = np.zeros((N, N), dtype=int)  # Initialize the nearest neighbors matrix

    # Create filenames for saving positions and velocities
    filenameX = f"./vicsek_cells/positionsN{N}L{L:.1f}r{r_0:.1f}v{v_0:.1f}eta{eta:.2f}.csv"
    filenameV = f"./vicsek_cells/velocitiesN{N}L{L:.1f}r{r_0:.1f}v{v_0:.1f}eta{eta:.2f}.csv"

    # Prepare data for saving
    X_data = np.array(X)
    V_data = np.array(V)

    # Save initial positions and velocities
    np.savetxt(filenameX, X_data, delimiter=",")
    np.savetxt(filenameV, V_data, delimiter=",")

    for k in range(nT):
        fly(X, Theta, V, dt, L, r_0, N, v_0, eta, cells, fixedRN)  # Evolve the system for one time step

        # Append current positions and velocities to data arrays
        X_data = np.array(X)
        V_data = np.array(V)

        # Append data to files
        with open(filenameX, "a") as fX, open(filenameV, "a") as fV:
            np.savetxt(fX, X_data, delimiter=",")
            np.savetxt(fV, V_data, delimiter=",")

        emptyCells(cells)  # Empty the cells for the next time step


# Generates a dataset of M realizations of the Vicsek model with N particles and nT time steps.
# Returns only the last state of each realization.
def DataSet2(N, M, dt, nT, r_0, v_0, eta, L, d):
    X_0 = np.random.rand(M, N, 2) * L  # Initialize random positions
    Theta_0 = np.random.rand(M, N) * 2 * np.pi  # Initialize random orientations
    V_0 = np.array([[[v_0 * np.cos(t), v_0 * np.sin(t)] for t in tt] for tt in Theta_0])  # Initialize velocities

    cells = CellGrid(d, L, r_0)  # Create the cell grid
    fixedRN = np.zeros((N, N), dtype=int)  # Initialize the nearest neighbors matrix
    np.fill_diagonal(fixedRN, 1)  # Each particle is a neighbor of itself

    # Evolve the first realization for nT time steps
    for _ in range(nT):
        fly(X_0[0], Theta_0[0], V_0[0], dt, L, r_0, N, v_0, eta, cells, fixedRN)
        emptyCells(cells)

    S = [momentConfig(X_0[0] / L, V_0[0] / v_0)]  # Calculate the moment configuration for the first realization

    # Evolve the remaining realizations for nT time steps
    for i in range(1, M):
        for _ in range(nT):
            fly(X_0[i], Theta_0[i], V_0[i], dt, L, r_0, N, v_0, eta, cells, fixedRN)
            emptyCells(cells)
        S.append(momentConfig(X_0[i] / L, V_0[i] / v_0))  # Calculate the moment configuration

    return X_0, V_0, S

# Returns a trajectory of length nT*M, sampled at every nT time steps
def dynamicDataset(N, M, dt, nT, r_0, v_0, eta, L, d):
    X = []  # Initialize an empty list to store positions
    X_0 = np.random.rand(N, 2) * L  # Initialize random positions
    Theta = []  # Initialize an empty list to store orientations
    Theta_0 = np.random.rand(N) * 2 * np.pi  # Initialize random orientations
    V = []  # Initialize an empty list to store velocities
    V_0 = np.array([[v_0 * np.cos(t), v_0 * np.sin(t)] for t in Theta_0])  # Initialize velocities

    cells = CellGrid(d, L, r_0)  # Create the cell grid
    fixedRN = np.zeros((N, N), dtype=int)  # Initialize the nearest neighbors matrix
    np.fill_diagonal(fixedRN, 1)  # Each particle is a neighbor of itself

    X.append(X_0.copy())  # Add the initial positions to the dataset
    Theta.append(Theta_0.copy())  # Add the initial orientations to the dataset
    V.append(V_0.copy())  # Add the initial velocities to the dataset

    # Evolve the system for nT time steps
    for _ in range(nT):
        fly(X[-1], Theta[-1], V[-1], dt, L, r_0, N, v_0, eta, cells, fixedRN)
        emptyCells(cells)

    S = [momentConfig(X_0 / L, V_0 / v_0)]  # Calculate the moment configuration for the initial state

    # Evolve the system for nT time steps for each realization
    for _ in range(1, M):
        X.append(X[-1].copy())  # Add the previous positions to the dataset
        Theta.append(Theta[-1].copy())  # Add the previous orientations to the dataset
        V.append(V[-1].copy())  # Add the previous velocities to the dataset
        for _ in range(nT):
            fly(X[-1], Theta[-1], V[-1], dt, L, r_0, N, v_0, eta, cells, fixedRN)
            emptyCells(cells)
        S.append(momentConfig(X[-1] / L, V[-1] / v_0))  # Calculate the moment configuration

    return X, V, S

def dynamicEnsemble(k, N, M, dt, nT, r_0, v_0, eta, L, d):
    ensemble = [dynamicDataset(N, M, dt, nT, r_0, v_0, eta, L, d)]
    for _ in range(1, k):
        ensemble.append(dynamicDataset(N, M, dt, nT, r_0, v_0, eta, L, d))
    return ensemble

# Calculated the optimum value of epsilon for a given number of particles, speed,
# interaction radius and noise intensity
def etaSet(N, M, dt, nT, r_0, v_0, etas, L, d, epsilons, D, maxoutdim=3):
    # Create filenames for saving data
    filename1 = f"./vicsekEvalsN{N}L{L:.1f}r{r_0:.2f}v{v_0:.2f}.csv"
    filename2 = f"./vicsekOrderN{N}L{L:.1f}r{r_0:.2f}v{v_0:.2f}.csv"
    filename3 = f"./vicsekSGEN{N}L{L:.1f}r{r_0:.2f}v{v_0:.2f}.csv"

    # Calculate and save data for each eta value
    for eta in etas:
        X, V, S = DataSet2(N, M, dt, nT, r_0, v_0, eta, L, d)
        for i in range(len(S)):
            for j in range(i):
                D[i, j] = np.linalg.norm(S[i] - S[j])
                D[j, i] = D[i, j]

        SGE, epsilonOpt = epsilonOptimum(D, epsilons)
        _, e_vals, _, _ = DiffMap(D, epsilonOpt, maxoutdim=maxoutdim)

        # Append data to files
        with open(filename1, "a") as f1, open(filename2, "a") as f2, open(filename3, "a") as f3:
            np.savetxt(f1, e_vals, delimiter=",")
            np.savetxt(f2, np.mean([np.sqrt(modes[1, i] ** 2 + modes[2, i] ** 2) for i in range(M)]), delimiter=",")
            np.savetxt(f3, SGE, delimiter=",")
            np.savetxt(f3, [epsilonOpt], delimiter=",")

# Calculates local minima, returns the largest value of epsilon that attains a minimum
# if it exists otherwise a fixed value
def epsilonOptimum(D, epsilons, epsilon_0=0.25):
    SGE = SGECriteria(D, epsilons)

    # Use scipy.signal.find_peaks to find local minima
    minima_indices, _ = find_peaks(-SGE)  # Find peaks of the inverted SGE
    if len(minima_indices) > 0:
        return SGE, epsilons[minima_indices[-1]]
    else:
        return SGE, epsilon_0