import  numpy as np
from sympy import Matrix, exp, sqrt

A = Matrix([[1,0,0],[0,1,1],[0,1,1]])
V1, D = A.diagonalize(normalize=True)

B = exp(A)
C = V1 * exp(D) * V1.T
a = 1/sqrt(2)
V = Matrix(
    [[1, 0, 0], 
    [0, a, -a], 
    [0, -a, a]]
    )

L = Matrix(
    [[1, 0, 0], 
    [0, 0, 0], 
    [0, 0, 2]]
    )

M1 = V * exp(L)
M2 = M1 * V.T




