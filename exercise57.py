import numpy as np
from quantumlib import X, I, I4, composite, diag, H, T, S, C, CNOT, Toffoli


M1 = composite(I, I, H)

M2 = CNOT(3, 1, 2)

M3 = composite(I, I, np.conj(T))

M4 = CNOT(3, 0, 2)

M5 = composite(I, I, T)

M6 = M2
M7 = M3
M8 = M4
M9 = composite(I, np.conj(T), T)
M10a = CNOT(3, 0, 1)
M10b = composite(I, I, H)
M11 = composite(I, np.conj(T), I)
M12 = M10a
M13 = composite(T, S, I)

M = M1.dot(M2).dot(M3).dot(M4).dot(M5).dot(M6).dot(M7).dot(M8).dot(M9).dot(M10a).dot(M10b).dot(M11).dot(M12).dot(M13)
T = Toffoli()
print(np.abs(T-M))