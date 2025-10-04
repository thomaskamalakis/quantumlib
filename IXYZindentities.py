import sympy as sp

"""
Define basic matrices
"""
I = sp.Matrix([[1,0],[0,1]])
X = sp.Matrix([[0, 1],[1, 0]])
Y = sp.Matrix([[0, -1j],[1j, 0]])
Z = sp.Matrix([[1, 0],[0, -1]])

H = 1 / sp.sqrt(2) * sp.Matrix([[1, 1], [1 ,-1]])
def f1(a):
    return (1+(-1)**a)/2

def f2(a):
    return (1-(-1)**a)/2

Xa = lambda a: f1(a) * I + f2(a) * X
Ya = lambda a: f1(a) * I + f2(a) * Y
Za = lambda a: f1(a) * I + f2(a) * Z
Ha = lambda a: f1(a) * I + f2(a) * H

a = sp.symbols(r'\alpha')

"""
Show that Y^a=Z^(1/2) X^a Z^(-1/2)
"""
M1 = sp.nsimplify(Ya(a) ,tolerance=1e-10, rational = True)
M2 = sp.nsimplify(Za(0.5) * Xa(a) * Za(-0.5),tolerance=1e-10, rational = True)
diff12 = sp.simplify(M1-M2)

"""
Show that X^a=H Z^a H
"""
M3 = H * Za(a) * H
M4 = Xa(a)
diff34 = sp.simplify(M3-M4)

"""
Show that H^a = Y^(1/4) Z^a Y^(-1/4)
"""

M5 = sp.simplify( Ya(1/4) * Za(a) * Ya(-1/4) )
M6 = sp.simplify(H**(a))


diff56 = sp.nsimplify(sp.simplify(M5-M6).evalf(), tolerance=1e-10)



