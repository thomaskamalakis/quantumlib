import numpy as np
from numbertheory.numlib import gcd

X = np.array([[0, 1],[1, 0]])
Y = np.array([[0, -1j],[1j, 0]])
Z = np.array([[1, 0],[0, -1]])

H = 1/np.sqrt(2) * np.array([[1,1],[1,-1]])
I = np.eye(2)
I4 = np.eye(4)
epi4 = np.exp(1j*np.pi/4)
T = np.array([[1, 0],[0, epi4]])
S = np.array([[1, 0],[0, 1j]])

def funm(M, func):
    """
    Returns the result of func over the diagonalizable matrix M
    """
    D, U = np.linalg.eig(M)
    D = func(D)
    Ut = U.conj().T
    return U @ np.diag(D) @ Ut
    
def expm(M):
    """
    Returns exp(M) where M is diagonizable
    """
    return funm(M, np.exp)

def rem_phase(M):
    """
    Sets the matrix element phase so that M[0,0] has zero phase
    """
    phi = np.angle(M[0,0])
    el = np.exp(-1j * phi)
    return M * el

def rotX(t):
    return np.array([
        [np.cos(t/2), -1j*np.sin(t/2)],
        [-1j*np.sin(t/2), np.cos(t/2)]
        ])

def rotY(t):
    return np.array([
        [np.cos(t/2), -np.sin(t/2)],
        [np.sin(t/2), np.cos(t/2)]
        ])

def rotZ(t):
    return np.array([
        [np.exp(-1j*t/2), 0],
        [0, np.exp(1j*t/2)]
        ])

def rotn(n, t):
    n = np.array(n)
    nx = n[0]
    ny = n[1]
    nz = n[2]
    return I * np.cos(t/2) - 1j * np.sin(t/2) * (nx * X + ny * Y + nz * Z)

def angles(U):
    u11 = U[0,0]
    u22 = U[1,1]
    u21 = U[1,0]
    p11 = np.angle(u11)
    p22 = np.angle(u22)
    p21 = np.angle(u21)
    
    c = 2 * np.arccos(np.abs(u11))
    a = (p11+p22)/2
    b = p21-p11
    d = p22-p21
    return a,b,c,d

def calcU(a,b,c,d):
    u11 =  np.cos(c/2) * np.exp(1j*(a-b/2-d/2))
    u12 = -np.sin(c/2) * np.exp(1j*(a-b/2+d/2))
    u21 =  np.sin(c/2) * np.exp(1j*(a+b/2-d/2))
    u22 =  np.cos(c/2) * np.exp(1j*(a+b/2+d/2))
    return np.array([[u11,u12],[u21,u22]])

            
def AXBXC(U):
    """
    Decompose a controlled C(U) using unitary and controlled-not gate
    """
    a, b, c, d = angles(U)
    A = np.matmul( rotZ(b), rotY(c/2) )
    B = np.matmul( rotY(-c/2), rotZ( -(d+b)/2) )
    C = rotZ( (d-b)/2 )
    return A,B,C
   
 
def diag(*args):
    sz = 0
    for M in args:
        sz += M.shape[0]
        
    D = np.zeros([sz,sz], dtype=complex)
    i = 0
    for M in args:
        sm = M.shape[0]
        D[i:i+sm,i:i+sm] = M
        i += sm
    
    return D
    
C = diag(I, X)  

def bin_array(num, m):
    """Convert a positive integer num into an m-bit bit vector"""
    return np.array(list(np.binary_repr(num).zfill(m))).astype(np.int8)

def bin_str_to_array(s):
    return np.array([int(d) for d in s ])

def bin_array_to_str(a):
    return ''.join([str(d) for d in a])

def to_dec(b):
    return int(''.join(map(str, b)), 2)
       
def CNOT(n, mc, mt):
    """
    Generalization of a CNOT gate with n input and output qubits
    The qubit indexed by mc is the control bit
    The qubit indexed by mt is the target it
    All other qubits go through unchanged
    Standard CNOT corresponds to n=2, mc=0, mt=1
    Returns the matrix for the CNOT gate operation
    """
    sz = 2 ** n
    C = np.zeros([sz,sz])
    
    for p in range(sz):
        bn = bin_array(p, n)
        if bn[mc] == 1:
            bn[mt] = 1 - bn[mt]
        q = to_dec(bn)
        C[p,q] = 1
    return C   


def basis_states(n):
    """
    Returns all basis states for n qubits
    """

    N = 2 ** n
    states = []
    for m in range(N):
        states.append( bin_array(m, n) )
    
    return states
        
def composite(*args):
    N = len(args)
    sz = 2 ** N
    D = np.zeros([sz,sz], dtype = complex)
    for n in range(sz):
        for m in range(sz):
            bn = bin_array(n, N)
            bm = bin_array(m, N)
            mul = 1            
            for i in range(N):
                mul *= args[i][ bn[i], bm[i] ]
            D[n, m] = mul 
    return D
 
def generalized_U(U, n, mt):
    """
    Return the matrix describing operating the 2x2 unitary matrix U on an array of n qubits 
    and where U only acts on the qubit indexed by mt
    """
    M = []
    I = np.eye(2)
    for i in range(n):
        if mt == i:
            M.append(U)
        else:
            M.append(I)
            
    return composite(*M)
             
def Toffoli(t = 2, N = 3):
    sz = 2 ** N
    T = np.zeros([sz,sz])
    for n in range(sz):
        bm = bin_array(n, N)
        if np.sum(bm) == N:
            bm[t] = 0
        elif ( np.sum(bm) == N-1 ) and (bm[t] == 0):
            bm[t] = 1
        m = to_dec(bm)
        T[n,m] = 1
    return T
    
def swap(a=1, b=2, N=3):
    sz = 2 ** N
    S = np.zeros([sz,sz])
    for n in range(sz):
        bm = bin_array(n, N)
        v = bm[a]
        bm[a] = bm[b]
        bm[b] = v
        m = to_dec(bm)
        S[n,m] = 1
    return S
    
def Friedkin():
    N = 3
    sz = 2 ** N
    F = np.zeros([sz,sz])
    for n in range(sz):
        bm = bin_array(n,N)
        if bm[0] == 1:
            v = bm[1]
            bm[1] = bm[2]
            bm[2] = v
        m = to_dec(bm)
        F[n,m] = 1
    return F
    
def level2unitary(U, m):
    """
    Estimate the unitary matrix that zeros the mth element of the first column of U
    """
    
    b = U[m,0]
    a = U[0,0]
    M = U.shape[0]
    
    if b == 0:
        return np.eye(M)
    
    Um = np.zeros([M,M], dtype = complex)
    ac = np.conj(a)
    bc = np.conj(b)
    n = np.sqrt( np.abs(a) ** 2 + np.abs(b) ** 2 )
    Um[0,0] = ac / n
    Um[0, m] = bc / n
    Um[m, 0] = b / n
    Um[m, m] = -a / n 
    for l in range(1, M):
        if l != m:
            Um[l,l] = 1
    
    return Um
            
def level2submatrices_1layer(U):
    """
    Find two-level unitary matrices of level M-1 where M is the dimension of U
    such that their product times U is equal to I
    """
    matrices = []
    M = U.shape[0]    
    U_tmp = np.copy(U)
    
    for m in range(1, M):
        Uc = level2unitary(U_tmp,m)
        U_tmp = Uc.dot(U_tmp)
        matrices.append(Uc)
        
    sub_matrix = U_tmp[1:,1:]
    
    return matrices, sub_matrix

def level2submatrices(U):
    """
    Decompose the a unitary matrix U in 2-level unitary sub-matrices
    """
    
    M = U.shape[0]    
    U_tmp = np.copy(U)
    matrices = []
    aug = 0
    
    for m in range(2, M):
        new_matrices, U_tmp = level2submatrices_1layer(U_tmp)
        print(U_tmp)
        if aug >=1 :
            I = np.eye(aug, dtype=complex)
            new_matrices = [diag(I, u) for u in new_matrices]
            
        matrices += new_matrices
        aug += 1
    
    I = np.eye(aug, dtype=complex)
    matrices.append(
        diag(I, U_tmp.conj().T)
        )
    
    U_appx = np.eye(M, dtype=complex)
    
    for v in matrices:
        U_appx = U_appx.dot(v.conj().T)
        
    return matrices, U_appx


def gray_count(initial,final):
       
    init_bin_array = bin_str_to_array(initial)
    final_bin_array = bin_str_to_array(final)
    xor = init_bin_array ^ final_bin_array
    elements = [ initial ]
    curr = np.copy(init_bin_array)
    
    for i, x in enumerate(xor):
        if x == 1:
            curr[i] = 1 - curr[i]
            s = bin_array_to_str(curr)
            elements.append( s )
    
    return elements

def qft(xm):
    """
    Perform quantum Fourier transform
    """
    N = xm.size
    return 1/np.sqrt(N) * np.fft.fft(xm)

def iqft(xm):
    """
    Perform inverse quantum Fourier transform
    """
    N = xm.size
    return np.sqrt(N) * np.fft.fft(xm)

def modU(xm, x, N):
    """
    Estimate the quantum state such that
    modU|y> = |xy mod N>    
    """
    
    zm = np.zeros(xm.size)
    
    if x >= N:
        raise ValueError('the second argument must be smaller than the third.')
        
    if gcd(x, N) != 1:
        raise ValueError('the second argument must be co-prime with the third.')
    
    for y in range(xm.size):
        
        i = np.mod(x * y, N, dtype=int)
        zm[i] = xm[y]
    
    return zm

class quantum_register:
    
    def __init__(self, sz = 1):
        self.combs = 2 ** sz
        self.sz = sz
        self.x = np.zeros( self.combs )
        indxs = np.arange(self.combs, dtype = int)        
        self.b = [
            np.binary_repr(el, width = sz) for el in indxs
            ]
        self.x[0] = 1
        
    
