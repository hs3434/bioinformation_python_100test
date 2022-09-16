with open("data/test.fq") as p:
    d = p.readlines()

i = 0
n = 0
for line in d:
    i += 1
    if i % 4 == 2:
        n += len(line) - 1
print(n)
