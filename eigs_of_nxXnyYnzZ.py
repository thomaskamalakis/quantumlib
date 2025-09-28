import sympy as s

nx, ny, nz = s.symbols('n_x, n_y, n_z', real=True)

X=s.Matrix([[0,1],[1,0]])
Y = s.Matrix([[0, -1j],[1j, 0]])
Z = s.Matrix([[1, 0],[0, -1]])

A = nx*X + ny*Y + nz*Z
P, D = A.diagonalize(normalize=True)
v1 = P[:,0]