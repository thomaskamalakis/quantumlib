import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc

rc('text.latex', preamble=r'\usepackage{amsmath}')
rc('text', usetex=True)
rc('font', size=20)
rc('lines',linewidth=1.5)
rc('savefig', dpi = 500)

phi = np.arccos( np.cos(np.pi/8) ** 2)

k = np.arange(0,10000)
phik = np.mod(k * phi, 2*np.pi)

plt.close('all')
fig = plt.hist(phik, bins=30, color='skyblue', edgecolor='black')
plt.xlabel(r'$\theta_k$')
plt.ylabel(r'$N(\theta_k)$')
plt.tight_layout()