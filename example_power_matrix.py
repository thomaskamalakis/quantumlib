import  numpy as np
from sympy import Matrix

A = Matrix([[1,0,0],[0,1,1],[0,1,1]])
V, D = A.diagonalize(normalize=True)

B = A ** 3
C = V * D **3 * V.TC



