import numpy as np

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

X = np.array([[0, 1],[1, 0]])
H = 1/np.sqrt(2) * np.array([[1,1],[1,-1]])
I = np.eye(2)
I4 = np.eye(4)
epi4 = np.exp(1j*np.pi/4)
T = np.array([[1, 0],[0, epi4]])
S = np.array([[1, 0],[0, 1j]])
            
def AXBXC(U):
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
        
    
    