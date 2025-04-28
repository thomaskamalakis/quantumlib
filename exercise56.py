import numpy as np
from quantumlib import X, I, AXBXC, angles

U = 1 / np.sqrt(2) * np.array([[1,1],[1,-1]])
D0 = np.array([[1, 0],[0, 1]])

X_numerical = np.conjugate(U).dot(D0).dot(U) 

D = np.array([[1, 0],[0, 1j]])

V =  np.conjugate(U).dot(D).dot(U)
X2 = V.dot(V)
V2 = 1j * (1-1j) / 2 * I + (1-1j) / 2 * X

A, B, C = AXBXC(V)