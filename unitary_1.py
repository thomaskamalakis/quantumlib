import numpy as np
from quantumlib import T, H, expm, rotX, rotY, X, Z, rotn, rem_phase

"""
Formula 4.74
"""

theta = 2 * np.arccos( np.cos(np.pi/8) ** 2)

nx = np.cos(np.pi/8) * np.sin(np.pi/8) / np.sin(theta/2) 
nz = nx
ny = np.sin(np.pi/8) ** 2 / np.sin(theta/2)

n = np.array([nx,ny,nz])

M1 = T @ H @ T @ H


M2 = expm(-1j*np.pi/8 * Z) @ expm(-1j*np.pi/8 * X)
M3 = rotn(n, theta)

print(rem_phase(M1))
print(rem_phase(M2))
print(rem_phase(M3))
