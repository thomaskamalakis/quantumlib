import numpy as np
from quantumlib import rotX, rotY, rotZ, angles, calcU, AXBXC, X

t = np.pi/6
RX = rotX(t)
a,b,c,d = angles(RX)
U = calcU(a,b,c,d)
A,B,C = AXBXC(RX)
U2 = A.dot(X).dot(B).dot(X).dot(C)

A2 = rotZ(-np.pi/2).dot( rotY(t/2) )
B2 = rotY(-t/2)
C2 = rotZ(np.pi/2)

R = rotZ(np.pi/2).dot(X).dot(rotY(t/2)).dot(X).dot(rotY(-t/2)).dot(rotZ(-np.pi/2))

R2 = X.dot(rotY(t/2)).dot(X).dot(rotY(-t/2))