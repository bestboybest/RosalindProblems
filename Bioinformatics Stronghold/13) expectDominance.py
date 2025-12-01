a = 19149
b = 19596
c = 16298
d = 16183
e = 17881
f = 16774

def expect(a, b, c, d, e, f):
    return 2* (a * 1 + b * 1 + c * 1 + d * 3/4 + e * 1/2 + f * 0)

print(expect(a, b, c, d, e, f))