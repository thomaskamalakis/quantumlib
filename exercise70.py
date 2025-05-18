import numpy as np

def myU(a,b, m = 0):
    ab = np.sqrt( np.abs(a) ** 2 +np.abs(b) ** 2 )
    ac = np.conj(a)
    bc = np.conj(b)
    U = np.zeros([3,3], dtype=complex)
    
    indices = []
    for n in range(3):
        if m!=n:
            indices.append(n)
    
    u = indices[0]
    d = indices[1]
    
    U[u,u] = ac/ab
    U[u,d] = bc/ab
    U[d,u] = b/ab
    U[d,d] = -b/ab
    U[m,m] = 1
    
    return U
 
U1 = myU(1+1j, 1-1j, m=0)
U2 = myU(1, 1, m=1)
U3 = myU(1, 1j, m=2)

U = U1.dot(U2).dot(U3)

print(U)
    
            
    