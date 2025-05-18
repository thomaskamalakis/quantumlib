import numpy as np
from quantumlib import level2submatrices

M = 0.5 * np.array([
        [1, 1, 1, 1],
        [1, 1j, -1, -1j],
        [1, -1, 1, -1],
        [1, -1j, -1, 1j]    
        ])

U, resU = level2submatrices(M) 

for i, u in enumerate(U):
    print('i=',i)
    print(u)
    print('\n')