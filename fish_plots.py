import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigs

#The next functions make a B-spline approximation of the trajectories, I guess these can be easily replaced with some python method, probably more efficient and cleaner

# Computes the B-spline approximation of a sequence of points X.
# X = [[x1,y1],[x2,y2],...] 
def BSplineAppex(X):
    n = len(X)  # Number of points

    c1 = [1.0 / 4]  # Initialize the first element of the c1 vector
    S = [(6 * X[1] - X[0]) / 4]  # Initialize the first element of the S vector (adjusting index for Python)

    # Recursively compute the remaining elements of S and c1
    for i in range(1, n - 3):  # Adjusting loop range for Python indexing
        S.append((6 * X[i + 1] - S[-1]) / (4 - c1[-1]))
        c1.append(1 / (4 - c1[-1]))

    # Compute the last element of S
    S.append((6 * X[n - 2] - X[n - 1] - S[-1]) / (4 - c1[-1]))  # Adjusting indices for Python
    c1.append(0)  # Append 0 to c1

    B = [S[-1]]  # Initialize the first element of the B vector

    # Recursively compute the remaining elements of B
    for i in range(len(S) - 1):
        B.append(S[-(i + 1)] - c1[-(i + 1)] * B[-1])

    # Return the reversed B vector
    return B[::-1]  # Using slicing for reversal

# Computes the lengths of consecutive sequences of integers in a vector I.
# Returns a vector of pairs [start_index, length] for each sequence.
def lenIs(I):
    if len(I) > 0:
        D = [[I[0] - 1, 1]]  # Initialize the first sequence (adjusting index for Python)
        for i in I[1:]:
            if i - (sum(D[-1])) == 1:  # Check if the current integer is consecutive to the previous sequence
                D[-1][1] += 1  # If consecutive, increase the length of the current sequence
            else:
                D.append([i - 1, 1])  # Otherwise, start a new sequence
        return D
    else:
        return [[0, 0]]  # Return [[0,0]] if I is empty

# Computes the cubic spline interpolation of a sequence of points X using BSplineAppex.
# Returns a vector of interpolated points.
def BSCubicInterpolator(X):
    B = BSplineAppex(X)  # Compute the B-spline approximation
    Q = [X[0]]  # Initialize the interpolated points with the first point in X (adjusting index for Python)
    n = len(X)  # Number of points

    # Add the first four interpolated points
    Q.append(X[0] + (B[0] - X[0]) / 3)  # Adjusting indices for Python
    Q.append(X[0] + 2 * (B[0] - X[0]) / 3)  # Adjusting indices for Python
    Q.append(X[1])  # Adjusting indices for Python
    Q.append(X[1])  # Adjusting indices for Python

    # Add the remaining interpolated points
    for i in range(1, n - 2):  # Adjusting loop range for Python indexing
        Q.append(B[i - 1] + (B[i] - B[i - 1]) / 3)
        Q.append(B[i - 1] + 2 * (B[i] - B[i - 1]) / 3)
        Q.append(X[i + 1])
        Q.append(X[i + 1])

    # Add the last four interpolated points
    Q.append(B[n - 3] + (X[n - 1] - B[n - 3]) / 3)  # Adjusting indices for Python
    Q.append(B[n - 3] + 2 * (X[n - 1] - B[n - 3]) / 3)  # Adjusting indices for Python
    Q.append(X[n - 1])  # Adjusting indices for Python

    return Q  # Return the interpolated points

# Interpolates a trajectory with missing data points.
# X = [[x1,y1],[x2,y2],...] where missing data points are replaced with [0,0].
# Is keeps track of the missing data points.
# The function fits a cubic spline to the non-zero data points and adds evenly spaced data points on the segments with missing data.
def trajectoryInterpolator(X):
    Is1 = np.where(np.linalg.norm(X, axis=1) == 0)[0]  # Find indices of missing data points using NumPy
    # Interpolate the non-zero data points using cubic spline
    X_nonzero = np.array([X[i] for i in range(len(X)) if i not in Is1])
    BS = BSCubicInterpolator(X_nonzero) 
    V = [3 * (BS[(i - 1) * 4 + 2] - BS[(i - 1) * 4 + 1]) for i in range(1, len(X_nonzero))]  # Calculate initial velocity estimates
    V.append(3 * (BS[-1] - BS[-2]))  # Add the last velocity estimate
    Is2 = lenIs(Is1)  # Compute the lengths of consecutive sequences of missing data points
    s = 0  # Initialize a counter for the number of missing data points processed

    # Iterate over the sequences of missing data points
    for i in range(len(Is2)):
        if Is2[i][0] == 0:  # If the missing data points are at the beginning
            for j in range(1, Is2[i][1] + 1):
                X[Is2[i][0] + (Is2[i][1] - j + 1)] = 2 * X[Is2[i][0] + (Is2[i][1] - j + 1) + 1] - X[Is2[i][0] + (Is2[i][1] - j + 1) + 2]  # Extrapolate backwards
                V.insert(0, V[0])  # Add a velocity estimate at the beginning
        elif sum(Is2[i]) == len(X):  # If the missing data points are at the end
            for j in range(1, Is2[i][1] + 1):
                X[Is2[i][0] + j] = 2 * X[Is2[i][0] + j - 1] - X[Is2[i][0] + j - 2]  # Extrapolate forwards
                V.append(V[-1])  # Add a velocity estimate at the end
        else:  # If the missing data points are in the middle
            h = 1.0 / (Is2[i][1] + 1)  # Calculate the spacing between interpolated points
            k = Is2[i][0] - s  # Calculate the index of the previous non-zero data point
            for j in range(1, Is2[i][1] + 1):
                # Interpolate the missing data points using cubic spline
                X[Is2[i][0] + j] = (
                    (1 - j * h) ** 3 * BS[(k - 1) * 4 + 1]
                    + 3 * (1 - j * h) ** 2 * j * h * BS[(k - 1) * 4 + 2]
                    + 3 * (1 - j * h) * (j * h) ** 2 * BS[(k - 1) * 4 + 3]
                    + (j * h) ** 3 * BS[(k - 1) * 4 + 4]
                )
                # Calculate the velocity estimate at the interpolated point
                V.insert(
                    Is2[i][0] + j,
                    3 * (1 - j * h) ** 2 * (BS[(k - 1) * 4 + 2] - BS[(k - 1) * 4 + 1])
                    + 6 * (1 - j * h) * j * h * (BS[(k - 1) * 4 + 3] - BS[(k - 1) * 4 + 2])
                    + 3 * (j * h) ** 2 * (BS[(k - 1) * 4 + 4] - BS[(k - 1) * 4 + 3]),
                )
        s += Is2[i][1]  # Update the counter for processed missing data points
    return np.array(V)  # Return the velocity estimates as a NumPy array

# Reads fish trajectory data from a file and interpolates missing data points.
# Returns the interpolated positions and velocities of the fish.
def fishPositions(N, m):
    filenameX = f"./fish/trajectories{N}_{m:02d}.csv"  # Construct the filename using f-string
    X = []  # Initialize an empty list to store the data
    with open(filenameX, "r") as Xio:  # Open the file for reading using 'with' statement
        for line in Xio:  # Iterate over each line in the file
            X.append([float(val) for val in line.strip().split(",")])  # Parse the line and append it to the list
    
    Y = [[X[(i - 1) * N + j] for i in range(1, len(X) // N + 1)] for j in range(N)]  # Reshape the data into a list of trajectories, adjusting indices for Python
    V = [trajectoryInterpolator(np.array(y)) for y in Y]  # Interpolate the trajectories, converting lists to NumPy arrays
    return [[Y[j][i] for j in range(N)] for i in range(len(Y[0]))], [[V[j][i] for j in range(N)] for i in range(len(V[0]))]  # Return the positions and velocities, adjusting indices for Python

# Renormalizes the data to the range [-1, 1] in both x and y dimensions.
def renormalisedData(X):
    xMin = np.min(X[:, 0])  # Find the minimum x value using NumPy
    yMin = np.min(X[:, 1])  # Find the minimum y value using NumPy
    xMax = np.max(X[:, 0])  # Find the maximum x value using NumPy
    yMax = np.max(X[:, 1])  # Find the maximum y value using NumPy
    return np.array([[(2 * x[0] - (xMin + xMax)) / (xMax - xMin), (2 * x[1] - (yMax + yMin)) / (yMax - yMin)] for x in X])  # Renormalize the data using NumPy

# Computes features for fish trajectories based on position and velocity.
def fishFeatures(N, m):
    X, _ = fishPositions(N, m)  # Get the fish positions
    #fishPositions gives us the velocity using spline interpolation, but it might be better to use a derivative with noise filtering
    V = noisyDerivative(X)  # Calculate the velocities
    X_2 = X[3:-3]  # Remove boundary points, adjusting indices for Python
    Xn = [renormalisedData(np.array(x)) for x in X_2]  # Renormalize the positions, converting lists to NumPy arrays
    Vn = [renormalisedData(np.array(x)) for x in V]  # Renormalize the velocities, converting lists to NumPy arrays
    MXV = [momentConfig(Xn[i], Vn[i]) for i in range(len(Xn))]  # Calculate the moment configurations
    return X_2, V, Xn, Vn, MXV  # Return the positions, velocities, renormalized positions and velocities, and moment configurations

# Computes features for fish trajectories based on position, velocity, and acceleration.
def fishAccelerationFeatures(N, m):
    X, _ = fishPositions(N, m)  # Get the fish positions
    V = noisyDerivative(X)  # Calculate the velocities
    A = noisyDerivative(V)  # Calculate the accelerations
    X_2 = X[6:-6]  # Remove boundary points, adjusting indices for Python
    V_2 = V[3:-3]  # Remove boundary points, adjusting indices for Python
    Xn = [renormalisedData(np.array(x)) for x in X_2]  # Renormalize the positions, converting lists to NumPy arrays
    Vn = [renormalisedData(np.array(v)) for v in V_2]  # Renormalize the velocities, converting lists to NumPy arrays
    An = [renormalisedData(np.array(a)) for a in A]  # Renormalize the accelerations, converting lists to NumPy arrays
    MXVA = [AccMomentConfig(Xn[i], Vn[i], An[i]) for i in range(len(Xn))]  # Calculate the moment configurations
    return X_2, V_2, A, Xn, Vn, An, MXVA  # Return the positions, velocities, accelerations, renormalized positions, velocities and accelerations, and moment configurations

# Computes z-score features for fish trajectories.
def fishZscoreFeatures(N, m, acc=False):
    X, _ = fishPositions(N, m)  # Get the fish positions
    V = noisyDerivative(X)  # Calculate the velocities
    X_1 = X[3:-3]  # Remove boundary points, adjusting indices for Python

    if acc:  # If acceleration features are requested
        A = noisyDerivative(V)  # Calculate the accelerations
        X_2 = X[6:-6]  # Remove boundary points, adjusting indices for Python
        V_2 = V[3:-3]  # Remove boundary points, adjusting indices for Python
        MXVA = [zScoreFeatures_acc(np.array(X_2[i]), np.array(V_2[i]), np.array(A[i])) for i in range(len(A))]  # Calculate the z-score features with acceleration, converting lists to NumPy arrays
        return X_2, V_2, A, MXVA  # Return the positions, velocities, accelerations, and z-score features
    else:  # If only position and velocity features are requested
        MXV = [zScoreFeatures(np.array(X_1[i]), np.array(V[i])) for i in range(len(V))]  # Calculate the z-score features, converting lists to NumPy arrays
        return X_1, V, MXV  # Return the positions, velocities, and z-score features

# Computes the noisy derivative of a sequence of points X using a central difference scheme.
def noisyDerivative(X):
    X = np.array(X)  # Convert X to a NumPy array for easier indexing
    return [(5 * (X[i + 1] - X[i - 1]) + 4 * (X[i + 2] - X[i - 2]) + X[i + 3] - X[i - 3]) / (2 ** 5) for i in range(3, len(X) - 3)]  # Calculate the derivative using a central difference scheme, adjusting indices for Python

#The rest of the functions are not needed, the plots can be reproduced using matplotlib once the diffusion map is done.