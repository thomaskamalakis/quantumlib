import sys
import os
import numpy as np

# Get the parent directory
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add parent directory to sys.path
sys.path.append(parent_dir)

import quantumlib as ql

X = 1+1j
Y = 2+5*1j
c_xy = np.sqrt(np.abs(X) ** 2 + np.abs(Y) ** 2)
X = X / c_xy
Y = Y / c_xy 
U = np.exp(1j)*np.array([[X, Y],[-np.conj(Y), np.conj(X)]])

phi_a = np.angle(U[0,0])
phi_b = np.angle(U[1,1])
phi = (phi_a+phi_b)/2

U2 = np.exp(-1j*phi) * U
x = U2[0 , 0]
y = U2[0 , 1]
L = np.abs(x)**2 / np.abs(y)**2
tan_phix = np.imag(x) / np.real(x)
tan_phiy = np.imag(y) / np.real(y)
tanq2 = np.sqrt(tan_phix ** 2 + (1+tan_phix**2) / L )
q = 2 * np.arctan(tanq2)
nz = -tan_phix / tanq2 

nx = 1
ny = 1/tan_phiy
nxy = np.sqrt(nx**2 + ny**2)
nx = nx / nxy * np.sqrt(1-nz**2)
ny = ny / nxy * np.sqrt(1-nz**2)

"""
nx and ny must be chosen so that the phase of -(ny+j*nx) corresponds to U2[1,0]
"""
if np.real(U2[1,0]) * ny < 0:
    ny = -ny
    nx = -nx
    
U3 = ql.rotn([nx, ny, nz], q) * np.exp(1j*phi)
