import sympy as sp
phi, theta, omega = sp.symbols(r'\phi, \theta, \omega', real=True)
"""
Συμβολικός υπολογισμός των διάφορων πινάκων περιστροφής στο τρισδιάστατο σύστημα συντεταγμένων
"""

"""
Βασικοί πίνακες περιστροφής γύρω από τους άξονες
"""
def Ry(p):
    cos = sp.cos(p)
    sin = sp.sin(p)
    return sp.Matrix([
        [cos, 0 ,  sin],
        [0,   1,    0],
        [-sin, 0,   cos]
        ])

def Rx(p):
    cos = sp.cos(p)
    sin = sp.sin(p)
    return sp.Matrix([
        [1, 0 ,  0],
        [0, cos, -sin],
        [0, sin,   cos]
        ])

def Rz(p):
    cos = sp.cos(p)
    sin = sp.sin(p)
    return sp.Matrix([
        [cos, -sin ,  0],
        [sin, cos, 0],
        [0, 0, 1]
        ])

"""
Υπολογισμός του πίνακα περιστροφής όπως το δείχνουμε στις διαφάνειες
"""
Rn = Rz(phi-sp.pi/2) * Rx(-theta) * Rz(omega) * Rx(theta) * Rz(sp.pi/2-phi) 
Rn = sp.simplify(Rn)

"""
Υπολογισμός από τον τύπο του Rodrigues
"""
nx = sp.sin(theta)*sp.cos(phi)
ny = sp.sin(theta)*sp.sin(phi)
nz = sp.cos(theta)

K = sp.Matrix([
    [0, -1 * nz, ny],
    [nz,   0,    -1*nx],
    [-1*ny,  nx, 0]
    ])

I3 = sp.eye(3)

Rn2 = I3 + sp.sin(omega) * K + (1-sp.cos(omega)) * K ** 2

"""
Διαφορά μεταξύ των δύο προσεγγίσεων
"""
diff2 = sp.simplify(Rn - Rn2)


