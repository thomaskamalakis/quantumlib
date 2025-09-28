import numpy as np
from quantumlib import array_to_latex

def my_norm(x,y):
    c = np.sqrt(np.abs(x)**2 + np.abs(y)**2)
    return x/c, y/c

M = 1/2 * np.array([
        [1, 1, 1, 1],
        [1, 1j, -1, -1j],
        [1, -1, 1, -1],
        [1, -1j, -1, 1j]    
        ])

A1 = np.array([
    [1/np.sqrt(2), 1/np.sqrt(2), 0 , 0],
    [1/np.sqrt(2), -1/np.sqrt(2), 0 , 0],
    [0, 0 , 1 , 0],
    [0, 0 , 0 , 1]])


C1  = A1.dot(M)
print( array_to_latex(C1) )

x, y = my_norm(1, C1[2,0] / C1[0,0])

A2 = np.array([
    [np.conj(x), 0, np.conj(y) , 0],
    [0, 1, 0 ,0 ],
    [-y, 0 , x , 0],
    [0, 0 ,0 , 1]])

C2 = A2.dot(C1)
print( array_to_latex(C2) )

x, y = my_norm(1, C2[3,0] / C2[0,0])
A3 = np.array([
    [np.conj(x), 0, 0 , np.conj(y)],
    [0, 1, 0 ,0 ],
    [0, 0 , 1, 0],
    [-y, 0 ,0 , x]])

C3 = A3.dot(C2)
print( array_to_latex(C3) )

C3s = C3[1:,1:]
x, y = my_norm(1, C3s[1,0] / C3s[0,0])
B3 = np.array([[np.conj(x), np.conj(y), 0],
                [-y, x, 0],
                [0, 0, 1]])

C4 = B3.dot(C3s)
print( array_to_latex(C4) )

x, y = my_norm(1, C4[2,0] / C4[0,0])
B5 = np.array([
    [np.conj(x), 0, np.conj(y)],
    [0, 1, 0],
    [-y, 0 , x]])

C5 = B5.dot(C4)
phase = np.angle(C5[0,0])
x, y = x * np.exp(1j*phase), y * np.exp(1j*phase) 

B5 = np.array([
    [np.conj(x), 0, np.conj(y)],
    [0, 1, 0],
    [-y, 0 , x]])
C5 = B5.dot(C4)


print( array_to_latex(C5) )

A4 = np.eye(4,4, dtype=complex)
A4[1:,1:] = B3

A5 = np.eye(4,4, dtype=complex)
A5[1:,1:] = B5

R=A5.dot(A4).dot(A3).dot(A2).dot(A1).dot(M)
print( array_to_latex(R) )



