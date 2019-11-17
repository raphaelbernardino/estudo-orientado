#!/usr/python3
from time import time
from primes_list import *

def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Newton's method
def isqrt(n):
    if n < 0:
        raise ValueError('Square root is not defined for negative numbers.')
    
    x = n
    y = (x + 1) // 2
    
    while y < x:
        x = y
        y = (x + n // x) // 2
    
    return x

def is_square(n):
    s = isqrt(n)
    return s * s == n


def factor(n, p=2):
    """
    Compute n-1 = 2^s * d for strong pseudoprime and strong lucas pseudoprime.
    """

    s = 0
    n -= 1

    while not n % p:
        s += 1
        n //= p

    return s, n


def jacobi(a, p):
    if (not p & 1) or (p < 0):
        raise ValueError('p must be a positive odd number.')

    if (a == 0) or (a == 1):
        return a

    a = a % p
    t = 1

    while a != 0:
        while not a & 1:
            a >>= 1
            if p & 7 in (3, 5):
                t = -t
        a, p = p, a
        if (a & 3 == 3) and (p & 3) == 3:
            t = -t
        a = a % p
    if p == 1:
        return t
    return 0


def selfridge(n):
    d = 5
    s = 1
    ds = d * s

    while True:
        if gcd(ds, n) > 1:
            return ds, 0, 0

        if jacobi(ds, n) == -1:
            return ds, 1, (1 - ds) // 4

        d += 2
        s *= -1
        ds = d * s


def lucas_sequence(n, u1, v1, u2, v2, d, q, m):
    k = q
    while m > 0:
        u2 = (u2 * v2) % n
        v2 = (v2 * v2 - 2 * q) % n
        q = (q * q) % n

        if m & 1:
            t1, t2 = u2 * v1, u1 * v2
            t3, t4 = v2 * v1, u2 * u1 * d
            u1, v1 = t1 + t2, t3 + t4

            if u1 & 1:
                u1 = u1 + n

            if v1 & 1:
                v1 = v1 + n

            u1, v1 = (u1 // 2) % n, (v1 // 2) % n
            k = (q * k) % n

        m = m >> 1

    return u1, v1, k


def strong_pseudoprime(n, base=2, s=None, d=None):
    if not s or not d:
        s, d = factor(n)

    x = pow(base, d, n)

    if x == 1:
        return True

    for i in range(s):
        if x == n - 1:
            return True
        x = pow(x, 2, n)
    return False


def strong_lucas_pseudoprime(n):
    d, p, q = selfridge(n)

    if p == 0:
        return n == d

    s, t = factor(n + 2)
    u, v, k = lucas_sequence(n, 1, p, 1, p, d, q, t >> 1)

    if (u == 0) or (v == 0):
        return True

    for i in range(1, s):
        v = (v * v - 2 * k) % n
        k = (k * k) % n
        if v == 0:
            return True

    return False


def baillie_psw(n, limit=100):
    if n == 2:
        return True

    if not n & 1:
        return False

    if n < 2 or is_square(n):
        return False

    bound = min(limit, isqrt(n))
    for i in range(3, bound, 2):
        if not n % i:
            return False

    return strong_pseudoprime(n, 2) and strong_lucas_pseudoprime(n)


def is_safe_prime(n):
    return is_prime(n) and is_prime((n - 1) // 2)


def is_prime(n):
    if n < 50:
        return n in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    return baillie_psw(n)


def test_primes():
    q = 0
    ini = time()
    for p in primes10000:
        if is_prime(p):
            q += 1
    end = time() - ini
    print('primes64: %d primes in %4.8f' % (q, end))

    q = 0
    ini = time()
    for p in primes64:
        if is_prime(p):
            q += 1
    end = time() - ini
    print('primes64: %d primes in %4.8f' % (q, end))

    q = 0
    ini = time()
    for p in primes128:
        if is_prime(p):
            q += 1
    end = time() - ini
    print('primes128: %d primes in %4.8f' % (q, end))

    q = 0
    ini = time()
    for p in primes256:
        if is_prime(p):
            q += 1
    end = time() - ini
    print('primes256: %d primes in %4.8f' % (q, end))

    q = 0
    ini = time()
    for p in primes512:
        if is_prime(p):
            q += 1
    end = time() - ini
    print('primes512: %d primes in %4.8f' % (q, end))


if __name__ == '__main__':
    test_primes()
