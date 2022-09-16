import pandas as pd

file1 = "data/test.fq"
file2 = "data/test.dup"
li = list()
n = 0
with open(file1) as p:
    while 1:
        d = p.readline()
        if d == "":
            break
        n += 1
        if n % 4 == 2:
            li.append(d)
b = pd.DataFrame(li)
c = b.drop_duplicates().index
n = 0
m = 0
d = list(["", "", "", ""])
with open(file1) as p, open(file2, mode="w") as p2:
    while 1:
        for v in range(4):
            d[v] = p.readline()
        if d[0] == "":
            break
        if n < len(c):
            if c[n] == m:
                p2.writelines(d)
                n += 1
        m += 1
