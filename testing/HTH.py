
import sympy as sp

phi = sp.symbols(r'\phi')
H = 1/sp.sqrt(2)*sp.Matrix([
    [1, 1],
    [1, -1]])

I = sp.eye(2)
T = sp.Matrix([[sp.exp(-1j*phi), 0], [0, sp.exp(1j*phi)]])

M = H * T * H

M = sp.simplify(M)

H1 = (I + H)/sp.sqrt(2)