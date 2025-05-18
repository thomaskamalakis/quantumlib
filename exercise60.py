from quantumlib import X, rotY
import numpy as np

M1 = rotY(np.pi/2)
M2 = X
M3 = rotY(-np.pi/2)
M = M1.dot(M2).dot(M3)
print(M)

RY = rotY(np.pi/4)
RYt = rotY(-np.pi/4)
M = RY.dot(X).dot(RY).dot(X).dot(RYt).dot(X).dot(RYt)
print(M)