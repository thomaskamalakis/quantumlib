from sympy import symbols, Matrix, sqrt, simplify, Eq, solve
a, b = symbols(r'\alpha \beta',real=True)
#b = sqrt(1-a**2)

X = Matrix([[0,1],[1,0]])
Z = Matrix([[1,0],[0,-1]])
I = Matrix([[1,0],[0,1]])

H = (b * X + a * Z)*a
P, D = H.diagonalize()
