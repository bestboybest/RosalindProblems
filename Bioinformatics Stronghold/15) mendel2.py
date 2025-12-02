# Aa x Bb -> 1/2 * 1/2 = 1/4 chance of AaBb progeny
# AA X Bb -> 1/2 * 1/2 = 1/4 chance of AaBb progeny
# Aa x BB -> same 1/4
# aa x Bb -> same 1/4

k = 7
n = 37

#4c1 * 1/4 * (1-1/4)^ + 4c2 * 1/16 + 4c3 * 1/64 + 4c4 * 1/64*4

def fact(a):
    if (a == 1 or a == 0):
        return 1
    return a * fact(a-1)

def mendel(k, n):
    #in K, there wil be 2**k progeny
    # atleast N out of 2**k have to be AaBb
    p = 2**k
    sum = 0
    for i in range(n, p+1):
        sum += (fact(p)/(fact(i) * fact(p-i))) * (1/(4 ** i)) * ((3/4)** (p-i))
    return sum

print(mendel(k, n))
