import sympy as sp

I = sp.Matrix([[1,0],[0,1]])

X = sp.Matrix([[0, 1],[1, 0]])
Y = sp.Matrix([[0, -1j],[1j, 0]])
Z = sp.Matrix([[1, 0],[0, -1]])

Rx = lambda o: I*sp.cos(o/2) - 1j * X * sp.sin(o/2)
Ry = lambda o: I*sp.cos(o/2) - 1j * Y * sp.sin(o/2)
Rz = lambda o: I*sp.cos(o/2) - 1j * Z * sp.sin(o/2)

def Rarb(o, m):
    mx = m[0]
    my = m[1]
    mz = m[2]
    
    A = (mx*X+my*Y+mz*Z)
    return I * sp.cos(o/2) -1j * A * sp.sin(o/2)

phi, theta = sp.symbols(r'\phi, \theta', real=True)
mx = sp.sin(theta) * sp.cos(phi)
my = sp.sin(theta) * sp.sin(phi)
mz = sp.cos(theta)

m = [mx, my, mz]

b,c,d = sp.symbols('b, c, d', real=True)
Rn = Rz(b) * Rarb(c, m) * Rz(d) 
Rn=sp.simplify(Rn)


