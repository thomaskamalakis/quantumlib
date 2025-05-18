import numpy as np
from quantumlib import basis_states

def single_out(e, indx):
    return np.array([el[indx] for el in e])

def components(e):
    n_components = e[0].size
    comps = []
    
    for i in range(n_components):
        comps.append(single_out(e, i))
        
    return tuple(comps)

def my_and(x,y,z):
    return x * y * z 

def my_and2(x,y,z):
    return x * y + x * z - x * (y ^ z)

def bar(n):
    return 1-n

def my_and_4v(x,y,z,w):
    return x * y * z * w
 
def my_and_4v2(x,y,z,w):
    return x*y + x*z + x*w - x*(y^z) - x*(y^w) - x*(z^w) + x*(y^z^w)

def my_and_5v(x,y,z,w,k):
    return x * y * z * w * k
 
def my_and_5v2(x,y,z,w, k):
    return ( 
       x*y + x*z + x*w + x*k 
      -x*(y^k) - x*(z^k) - x*(w^k)
      -x*(y^z) - x*(y^w) 
      -x*(z^w)
      +x*(y^k^z) + x*(y^k^w) + x *(z^w^k) + x*(y^z^w)
      -x*(y^z^w^k)      
     )

def print_truth_table( bs, v ):
    for i, b in enumerate(bs):
        for j in range(b.size):
            print(b[j], end="")
        print(' :', end="")
        print('%3d' % v[i])

    
xyzwk = basis_states(5)
x,y,z,w,k = components(xyzwk)

f1 = my_and_5v(x,y,z,w,k)
f2 = my_and_5v2(x,y,z,w,k)
print(f1)
print(f2)
print_truth_table(xyzwk, f2)