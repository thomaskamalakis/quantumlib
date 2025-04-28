from sympy import Function, Matrix, eye, vector, symbols, sin, cos, simplify, I, sqrt

nx, ny, nz = symbols('n_x n_y n_z', real=True)
kx, ky, kz = symbols('k_x k_y k_z', real=True)
theta = symbols(r'\theta',real=True)

#n = Matrix([[nx],[ny],[nz]])
n = Matrix([[nx],[ny],[sqrt(1-nx**2-ny**2)]])
k = Matrix([[kx],[ky],[kz]])
vrot = k*cos(theta) + n.cross(k)*sin(theta) + n*(n.dot(k))*(1-cos(theta))

X = Matrix([[0, 1],
            [1, 0]])

Y = Matrix([[0, -I],
            [I, 0]])

Z = Matrix([[1, 0],
            [0, -1]])

Mn = n[0]*X + n[1]*Y + n[2]*Z
Mk = kx*X + ky*Y + kz*Z

Mnk = cos(theta/2)**2 * Mk + I*cos(theta/2)*sin(theta/2) * Mk * Mn - I*cos(theta/2)*sin(theta/2) * Mn * Mk + sin(theta/2)**2 * Mn * Mk * Mn
Mnk = simplify(Mnk)

vx = vrot[0]
vy = vrot[1]
vz = vrot[2]
Mv = vx*X + vy*Y + vz*Z

M = Mnk - Mv
M0 = simplify(M)


