from numlib import composite, gcd, CRM_solve
import numpy as np

mj = [21,11]
xj = [1, -1]
M = np.prod(mj)

"""
Find a non-trivial solution of x^2 = 1(mod N)
"""
print('M = %d' %M)
x = CRM_solve(xj, mj)
print('x = %d is a solution' %x)
a2 = np.mod(x ** 2, M)
print('x^2 = %d mod %d = %d' %(x ** 2, M, a2) )

a = np.mod(x, M)
print('a = %d' %a)
g1 = gcd(M, a+1)
g2 = gcd(M, a-1)
print('g1 and g2 are equal to %d and %d respectively.' %(g1, g2))

