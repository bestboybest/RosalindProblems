def mortal(n, m):
    B = [0] * (n + 1)
    B[1] = 1

    for i in range(2, n + 1):
        for k in range(2, m + 1):
            if i - k > 0:
                B[i] += B[i - k]

    total = 0
    for k in range(m):
        if n - k > 0:
            total += B[n - k]

    return total


print(mortal(90, 19))
