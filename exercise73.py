import numpy as np
import matplotlib.pyplot as plt
from quantumlib import S, X, rem_phase, Z

G = 1/4*(3*S + X @ S @ X)
N = G.conj().T @ G
G0 = rem_phase(G)

print(G0)


G1 = 1/4 * (S + X @ S @ X )
G1_0 = rem_phase(G1)
print(G1_0)

G2 = 1/4 * (S - X @ S @ X )
G2_0 = rem_phase(G2)
print(G2_0)
