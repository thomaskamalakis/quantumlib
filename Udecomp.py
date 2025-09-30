import sys
import os
import numpy as np
import quantumlib as ql

X = 1+1j
Y = 2+5*1j
c_xy = np.sqrt(np.abs(X) ** 2 + np.abs(Y) ** 2)
X = X / c_xy
Y = Y / c_xy 
U = np.array([[X, Y],[-np.conj(Y), np.conj(X)]])

phi, q, [nx, ny, nz] = ql.decompU(U)
U3 = ql.rotn([nx, ny, nz], q) * np.exp(1j*phi)
