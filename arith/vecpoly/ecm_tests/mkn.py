import random, sys


def is_prime(n):
    if n < 2:
        return False
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def nextprime(n):
    n |= 1
    while not is_prime(n):
        n += 2
    return n


# usage: mkn.py seed pbits nbits
seed, pbits, nbits = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
random.seed(seed)
p = nextprime(random.getrandbits(pbits) | (1 << (pbits - 1)))
q = nextprime(random.getrandbits(nbits - pbits) | (1 << (nbits - pbits - 1)))
print(p, p * q)
