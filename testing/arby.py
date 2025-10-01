import sys
import os
import numpy as np

# Get the parent directory
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add parent directory to sys.path
sys.path.append(parent_dir)

import quantumlib as ql

phi=np.pi/6
theta=np.pi/3


nx = np.sin(theta)*np.cos(phi)
ny = np.sin(theta)*np.sin(phi)
nz = np.cos(theta)

n = ql.spherical_to_cartesian(1, theta, phi)
(nx, ny, nz) = n

n = np.array([[nx],[ny],[nz]])
nrot1 = ql.Rx( theta-np.pi/2 ) @ ql.Rz(np.pi/2 - phi) @ n 

"""
Construct a perpendicular vector to n
"""
e = np.array([1,0,0])
m = np.cross(e,[nx,ny,nz])
m = m / np.linalg.norm(m)

mrot1 = ql.Rx( theta-np.pi/2 ) @ ql.Rz(np.pi/2 - phi) @ m.T
theta_m = np.arccos(mrot1[2])

mrot2 = ql.Ry(-theta_m) @ mrot1

