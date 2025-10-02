import sympy as sp

"""
Define basic matrices
"""
I = sp.Matrix([[1,0],[0,1]])
X = sp.Matrix([[0, 1],[1, 0]])
Y = sp.Matrix([[0, -1j],[1j, 0]])
Z = sp.Matrix([[1, 0],[0, -1]])
H = 1/sp.sqrt(2) * sp.Matrix([[1, 1],[1, -1]])

"""
Rotation matrices
"""
Rx = lambda o: I*sp.cos(o/2) - 1j * X * sp.sin(o/2)
Ry = lambda o: I*sp.cos(o/2) - 1j * Y * sp.sin(o/2)
Rz = lambda o: I*sp.cos(o/2) - 1j * Z * sp.sin(o/2)

"""
Check properties of Hadamard matrices
Define the square root of H (H2) and the inverse square toot of H (H1)
"""

H2 = (1+1j)/2 * I + (1-1j)/2 * H
H1 = (1-1j)/2 * I + (1+1j)/2 * H

"""
Components of a unitary vector
"""
nx, ny, nz = sp.symbols(r'n_x, n_y, n_z')
M2 = sp.simplify(H2 * (nx * X + ny * Y + nz * Z) * H1)


mx =(nx - ny * sp.sqrt(2) + nz)/2
my = sp.sqrt(2)/2 * (nx-nz)
mz = (nx + ny * sp.sqrt(2) +nz)/2
M3 = mx * X + my * Y + mz *Z
diff = sp.simplify(M2-M3)

diffX = sp.simplify(H2 * X * H1 -0.5 * (X + sp.sqrt(2)*Y + Z))
diffY = sp.simplify(H2 * Y * H1 -0.5 * sp.sqrt(2) * (Z-X))
diffZ = sp.simplify(H2 * Z * H1 -0.5 * (X - sp.sqrt(2)*Y + Z))


mn = sp.simplify(mx * nx + my * ny + mz * nz)

