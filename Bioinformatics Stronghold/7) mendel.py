# k homozygous dominant, m heterozygous , n homozygous recessive
# We want probability that individual shows dominant phenotype

# Probability = 1 - (probability of recessive phenotype)

k = 26
m = 24
n = 17

def dominant(k, m, n):
    Recessive = m * (m-1)/2 * 1/4 + m * n * 1/2 + n * (n-1)/2 * 1
    Total = (k + m + n) * (k + m + n - 1)/2
    pDominant = 1 - (Recessive/Total)
    return pDominant

print(dominant(k, m, n))

