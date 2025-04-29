from quantumlib import swap, Toffoli, Friedkin, CNOT

M1 = CNOT(3,1,2)
M2 = swap()
M3 = Toffoli()
M4 = swap()
M5 = CNOT(3,1,2)

M = M1.dot(M2).dot(M3).dot(M4).dot(M5)
F = Friedkin()

print(M-F)