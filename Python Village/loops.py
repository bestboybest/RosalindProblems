a = 4144
b = 8590
sum = 0

for i in range(a, b + 1):
    if (i % 2 != 0):
        sum += i

print(sum)