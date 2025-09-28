import  numpy as np

A = np.matrix([[1,0,0],
               [0,1,1],
               [0,1,1]])
B = np.matrix([[2,1,1],
               [1,1,0],
               [0,0,1]])

C = A.dot(B)

D,U = np.linalg.eig(C)
D = np.diag(D)

res1 = C.dot(U)-U.dot(D)
print(res1)




