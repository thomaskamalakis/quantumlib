from quantumlib import AXBXC, X, angles
import numpy as np

fk = np.pi/6

U = np.array([[1, 0],
              [0, np.exp(1j*fk)]])
a,b,c,d = angles(U)

A, B, C = AXBXC(U)

U1 = np.exp(1j*a) * A.dot(X).dot(B).dot(X).dot(C)

