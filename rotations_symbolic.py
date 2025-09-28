import sympy as sp
"""
Συμβολικός υπολογισμός των διάφορων πινάκων περιστροφής σε επίπεδο qubit
"""

phi, theta, omega = sp.symbols(r'\phi, \theta, \omega', real=True)

"""
Ορισμός βασικών πινάκων
"""
I = sp.Matrix([[1,0],[0,1]])

X = sp.Matrix([[0, 1],[1, 0]])
Y = sp.Matrix([[0, -1j],[1j, 0]])
Z = sp.Matrix([[1, 0],[0, -1]])

Rx = lambda o: I*sp.cos(o/2) - 1j * X * sp.sin(o/2)
Ry = lambda o: I*sp.cos(o/2) - 1j * Y * sp.sin(o/2)
Rz = lambda o: I*sp.cos(o/2) - 1j * Z * sp.sin(o/2)

"""
Υπολογισμός του πίνακα περιστροφής όπως το δείχνουμε στις διαφάνειες
"""
Rn = Rz(phi-sp.pi/2) * Rx(-theta) * Rz(omega) * Rx(theta) * Rz(sp.pi/2-phi) 
Rn = sp.simplify(Rn)

"""
Υπολογισμός του πίνακα περιστροφής όπως υπολογίζεται από το θεωρητικό αποτέλεσμα
"""
nx = sp.sin(theta)*sp.cos(phi)
ny = sp.sin(theta)*sp.sin(phi)
nz = sp.cos(theta)


A = (nx*X+ny*Y+nz*Z)
A = sp.simplify(A) 
Rn2 = I * sp.cos(omega/2) - 1j * A * sp.sin(omega/2)

"""
Διαφορά μεταξύ των δύο προσεγγίσεων
"""
diff = sp.simplify(Rn-Rn2)



