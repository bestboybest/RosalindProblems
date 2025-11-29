f = open('./Input/rosalind_ini5.txt', 'r')

i = 1
for line in f:
    if (i % 2 == 0):
        print(line.strip())
    i += 1

