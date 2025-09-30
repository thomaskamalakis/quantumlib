import sympy as sp
mx, my, mz = sp.symbols('m_x, m_y, m_z', real=True)
b, c, d = sp.symbols('b, c, d', real=True)


I = sp.Matrix([[1,0],[0,1]])

X = sp.Matrix([[0, 1],[1, 0]])
Y = sp.Matrix([[0, -1j],[1j, 0]])
Z = sp.Matrix([[1, 0],[0, -1]])

Rx = lambda o: I*sp.cos(o/2) - 1j * X * sp.sin(o/2)
Ry = lambda o: I*sp.cos(o/2) - 1j * Y * sp.sin(o/2)
Rz = lambda o: I*sp.cos(o/2) - 1j * Z * sp.sin(o/2)
Rm = lambda o: I*sp.cos(o/2) - 1j * (mx*X+my*Y+mz*Z) * sp.sin(o/2)

R = Rz(b) * Rm(c) * Rz(d)
R = sp.simplify(R)