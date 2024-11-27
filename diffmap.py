import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigs

# Computes the mean and covariance matrix of concatenated state (X) and velocity (U) vectors.
# Returns a flattened vector of the means and covariance matrix elements.
# X = [[x1,y1],[x2,y2],...], U=[[u1,v1],[u2,v2],...]
def momentConfig(X, U):
    Z = np.concatenate((X, U), axis=1)  # Concatenate state and control vectors
    meanZ = np.mean(Z, axis=0)  # Calculate the mean of the concatenated vectors
    covZ = np.cov(Z, rowvar=False)  # Calculate the covariance matrix of the concatenated vectors
    # Return a vector of means and covariance matrix elements
    return np.array([
        meanZ[0], meanZ[1], meanZ[2], meanZ[3], covZ[0, 0], covZ[0, 1], covZ[0, 2], covZ[0, 3], covZ[1, 1],
        covZ[1, 2], covZ[1, 3], covZ[2, 2], covZ[2, 3], covZ[3, 3]
    ])

# Computes the mean and covariance matrix of concatenated state (X), velocity (U), and acceleration (A) vectors.
# Returns a flattened vector of the means and covariance matrix elements.
# X = [[x1,y1],[x2,y2],...], U=[[u1,v1],[u2,v2],...], A=[[ax1,ay1],[ax2,ay2],...]
def AccMomentConfig(X, U, A):
    Z = np.concatenate((X, U, A), axis=1)  # Concatenate state, control, and acceleration vectors
    meanZ = np.mean(Z, axis=0)  # Calculate the mean of the concatenated vectors
    covZ = np.cov(Z, rowvar=False)  # Calculate the covariance matrix of the concatenated vectors
    # Return a vector of means and covariance matrix elements
    return np.array([
        meanZ[0], meanZ[1], meanZ[2], meanZ[3], meanZ[4], meanZ[5], covZ[0, 0], covZ[0, 1], covZ[0, 2], covZ[0, 3],
        covZ[0, 4], covZ[0, 5], covZ[1, 1], covZ[1, 2], covZ[1, 3], covZ[1, 4], covZ[1, 5], covZ[2, 2], covZ[2, 3],
        covZ[2, 4], covZ[2, 5], covZ[3, 3], covZ[3, 4], covZ[3, 5], covZ[4, 4], covZ[4, 5], covZ[5, 5]
    ])

# Computes the z-scores of a given vector X.
def zScores(X):
    meanX = np.mean(X, axis=0)  # Calculate the mean of X
    stdX = np.std(X, axis=0)  # Calculate the standard deviation of X
    return (X - meanX) / stdX  # Calculate and return the z-scores

# Computes z-score features for state (X) and control (U) vectors.
def zScoreFeatures(X, U):
    zX = zScores(X)  # Calculate z-scores for X
    zU = zScores(U)  # Calculate z-scores for U
    return momentConfig(zX, zU)  # Return the moment configuration of the z-scored vectors

# Computes z-score features for state (X), control (U), and acceleration (A) vectors.
def zScoreFeatures_acc(X, U, A):  # Renamed to avoid conflict with previous function
    zX = zScores(X)  # Calculate z-scores for X
    zU = zScores(U)  # Calculate z-scores for U
    zA = zScores(A)  # Calculate z-scores for A
    return AccMomentConfig(zX, zU, zA)  # Return the moment configuration of the z-scored vectors

# Computes the Diffusion Map matrix L using the distance matrix D.
def DiffMap(D, epsilon, maxoutdim=2, maxiter=1000):
    Dmax = np.max(D)  # Maximum value in D
    K = np.exp(-(D**2) / (Dmax**2 * epsilon))  # Compute the kernel matrix

    k = np.diag(np.sum(K, axis=1))  # Compute the diagonal matrix of row sums
    L = np.linalg.inv(k) @ K  # Compute the normalized Laplacian matrix

    # Use eigh for symmetric matrices (should be faster and more stable)
    eigenvalues, e = eigh(L, eigvals=(len(L) - maxoutdim, len(L) - 1))  
    eigenvalues = eigenvalues[::-1]  # Reverse eigenvalues to descending order
    e = e[:, ::-1]  # Reverse eigenvectors accordingly
    Y = (eigenvalues * e).T  # Compute the diffusion map embedding

    return L, eigenvalues, e, Y  # Return Laplacian, eigenvalues, eigenvectors and embedding

# Computes the distance matrix from a set of features.
def distance_matrix(features):
    n = len(features)  # Number of features
    D = np.zeros((n, n))  # Initialize distance matrix
    for i in range(n):
        for j in range(i):  # Python uses 0-based indexing
            D[i, j] = np.linalg.norm(features[i] - features[j])  # Calculate Euclidean distance
            D[j, i] = D[i, j]  # Ensure symmetry
    return D  # Return distance matrix

# Computes the Semigroup Error (SGE) criterion for a given distance matrix D and epsilon.
def SGECriterion(D, epsilon):
    Dmax = np.max(D)  # Maximum distance

    K1 = np.exp(-(D**2) / (Dmax**2 * epsilon))  # Compute kernel matrix with epsilon
    k1 = np.diag(np.sum(K1, axis=1))  # Compute diagonal matrix of row sums
    L1 = np.linalg.inv(k1) @ K1  # Compute normalized Laplacian

    K2 = np.exp(-(D**2) / (Dmax**2 * (2 * epsilon)))  # Compute kernel matrix with 2*epsilon
    k2 = np.diag(np.sum(K2, axis=1))  # Compute diagonal matrix of row sums
    L2 = np.linalg.inv(k2) @ K2  # Compute normalized Laplacian

    return np.sqrt(np.trace((L1 @ L1 - L2) @ (L1 @ L1 - L2).T))  # Return SGE criterion using trace


# Computes the Semigroup Error (SGE) criterion for a given distance matrix D and epsilon.
def SGECriterion2(D, epsilon):
    Dmax = np.max(D)  # Maximum distance

    K1 = np.exp(-(D**2) / (Dmax**2 * epsilon))  # Compute kernel matrix with epsilon
    k1 = np.diag(np.sum(K1, axis=1))  # Compute diagonal matrix of row sums
    L1 = np.linalg.inv(k1) @ K1  # Compute normalized Laplacian

    K2 = np.exp(-(D**2) / (Dmax**2 * (2 * epsilon)))  # Compute kernel matrix with 2*epsilon
    k2 = np.diag(np.sum(K2, axis=1))  # Compute diagonal matrix of row sums
    L2 = np.linalg.inv(k2) @ K2  # Compute normalized Laplacian

    # Using eigh for symmetric matrices
    w, _ = eigh((L1 @ L1 - L2) @ (L1 @ L1 - L2).T, eigvals=(len(L1)-1, len(L1)-1))  
    return np.sqrt(w[0])  # Return SGE criterion using eigenvalues

# Computes the normalized SGE criteria for a range of epsilon values.
def SGECriteria(D, epsilons):
    SGE = np.array([SGECriterion(D, epsilon) for epsilon in epsilons])  # Compute SGE for each epsilon
    maxSGE = np.max(SGE)  # Find maximum SGE value
    return SGE / maxSGE  # Normalize and return SGE values

# Computes the normalized SGE criteria using eigenvalues for a range of epsilon values.
def SGECriteria2(D, epsilons):
    SGE = np.array([SGECriterion2(D, epsilon) for epsilon in epsilons])  # Compute SGE for each epsilon
    maxSGE = np.max(SGE)  # Find maximum SGE value
    return SGE / maxSGE  # Normalize and return SGE values

# Computes the connectivity criterion for a given distance matrix D and epsilon.
def ConnectivityCriterion(D, epsilon):
    Dmax = np.max(D)  # Maximum distance
    return np.sum(D / Dmax < epsilon) / len(D.flatten())  # Return connectivity criterion

# Computes the connectivity criteria for a range of epsilon values.
def ConnectivityCriteria(D, epsilons):
    return np.array([ConnectivityCriterion(D, epsilon) for epsilon in epsilons])  # Return connectivity criteria for each epsilon