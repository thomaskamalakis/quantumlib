from numlib import order_mod, gcd
import numpy as np

M = 11
a = 21

r = order_mod(a, M)
print('order of %d modulo %d is %d' %(a, M, r))
print('%d^%d = %d = %d modulo %d' %(a,r,a**r,np.mod(a**r,M),M))

