n = 31
k = 2

#Recurrence relation: Fn = Fn-1 + kFn-2

def fibo(n, k):
    if (n == 1):
        return 1
    if (n == 2):
        return 1
    return fibo(n-1, k) + k * fibo(n-2, k)

print(fibo(n, k))