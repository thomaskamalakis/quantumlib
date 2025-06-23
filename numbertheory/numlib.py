import numpy as np

def gcd(a,b):
    if a < b:
        a , b = b, a

    r1 = a
    r2 = b    

    while True:

        k = np.floor(r1/r2)
        r3 = r1 - k *r2
    
        if r3 == 0:
            break
        
        r1 = r2
        r2 = r3
    
    return r2

def are_co_prime(a,b):
    return gcd(a,b) == 1

def min_comb(a,b):
    if a < b:
        a , b = b, a
        reverse = True
    else:
        reverse = False

    r1 = a
    r2 = b    
    rr = [a, b]
    kk = []
    while True:

        k = np.floor(r1/r2)
        r3 = r1 - k *r2
        
        rr.append(r3)
        kk.append(k)
        
        if r3 == 0:
            break
        
        r1 = r2
        r2 = r3
        
    rr = np.array(rr)
    kk = np.array(kk)
    N = kk.size
    if N==1:
        if reverse:
            return 1, 0, b
        else:
            return 0, 1, b
        
    gcd = rr[N]
    A = 1
    B = -kk[N-2]
    for m in range(N-3, -1, -1):
       A, B = B, A - kk[m] * B
       
    if reverse:
        B, A = A, B
    return A, B, gcd

def inverse_modulo(a, n):
    A, B, gcd = min_comb(a, n)
    if gcd == 1:
        return A
    else:
        raise ValueError('%d does not have an inverse (mod) %d (gcd=%d)' %(a,n,gcd))


def prime_factors(n):
    factors = []

    # Strip out any prime factors that are 2.
    if n >= 2:
        while n % 2 == 0:
            factors.append(2)
            n //= 2
  
    # Check odd factors from 3 up to the square root of n.
    i = 3
    while i < n // i:
        while n % i == 0:
            factors.append(i)
            n //= i
        i += 2

    # If n is still greater than 1, it must be a prime factor.
    if n > 1:
        factors.append(n)
    
    return factors

def CRM_solve(xj, mj):
    """
    Chinese remainder theorem
    """
    xj = np.array(xj)
    mj = np.array(mj)
    
    N = mj.size
    for p in range(N):
        for q in range(N):
            if p!=q:
                g = gcd(mj[p], mj[q])
                if g != 1:
                    raise ValueError('the modulo parameters must all be co-prime. This is not the case for %d %d which have GCD equal to %d' %(mj[p], mj[q], g))
        
    Mj = np.zeros(N)
    Nj = np.zeros(N)
    M = np.prod(mj)
    
    for p in range(N):
        Mj[p] = M/mj[p]
        Nj[p] = inverse_modulo(Mj[p], mj[p])
        
    return np.sum(xj * Mj * Nj)

def CRM_RHS(x, mj):
    """
    Right hand side of the Chinese remainder theorem
    """
    xj = np.zeros(mj.size)
    for p in range(mj.size):
        xj[p] = np.mod(x, mj[p])
    
    return xj

def composite(pi,ai):
    """
    Estimation of a number from its primary number decomposition
    """
    pi = np.array(pi)
    ai = np.array(ai)
    return np.prod(pi**ai)

def order_mod(a, N):
    
    r = 1
    ar = a
    if np.gcd(a,N) != 1:
        raise ValueError('%d is not co-prime with %d, greatest commmon divisor is: %d' %(a,N,gcd(a,N)))
    
    while np.mod(ar, N) != 1:
        r += 1
        ar = ar * a
        
    return r

def is_even(n):
    return np.mod(n, 2) == 0

def is_odd(n):
    return np.mod(n, 2) == 1

def divides(a, b):
    return np.mod(b, a) == 0

def is_modulo(a, b, N):
    return np.mod(a, N) == np.mod(b, N)

def find_factor(N, select=None):
    while True:
        if select:
            y = select
        else:
            y = np.random.randint(N)
        g = np.gcd(y, N)
        if g == 1:
            r = order_mod(y, N)
            if is_even(r):
                yr = y ** int(r/2)
                if (not is_modulo(yr, -1, N)) and (not is_modulo(yr, 1, N)):
                    gg = [np.gcd(yr-1, N), np.gcd(yr+1,N)]
                    return gg

        
                
                