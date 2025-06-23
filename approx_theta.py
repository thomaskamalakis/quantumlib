import numpy as np
import matplotlib.pyplot as plt

Ns = np.arange(10, 10000, 1)

ns = []
theta = 2 * np.arccos( np.cos(np.pi/8) ** 2)
errors = []

for N in Ns:
    d = 2 * np.pi / N
    n = 1
    phi = np.mod( n * theta, 2*np.pi)
    error = np.abs(phi)
    
    while error > d:
        phi = np.mod( n * theta, 2*np.pi)
        error = np.abs(phi)
        n += 1
        
    ns.append(n)
    errors.append(error)

plt.close('all')
plt.figure()
plt.plot(Ns, ns)

plt.figure()
plt.plot(ns, errors)




