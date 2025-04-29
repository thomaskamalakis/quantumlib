from quantumlib import swap, Toffoli, Friedkin

M1 = Toffoli()
M2 = swap()
M3 = Toffoli()
M4 = swap()
M5 = Toffoli()

M = M1.dot(M2).dot(M3).dot(M4).dot(M5)
F = Friedkin()

print(M-F)