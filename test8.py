import numpy as np
import pandas as pd

file1 = "data/test.fq"
file2 = "data/test.cg"
a = np.empty((0, 150))
n = 0
with open(file1) as p:
    while 1:
        d = p.readline()
        if d == "":
            break
        n += 1
        if n % 4 == 2:
            tmp = np.frombuffer(str(d[0:-1]).encode("utf-8"), dtype='S1')
            a = np.row_stack((a, tmp))

b = pd.DataFrame(a)
c = pd.DataFrame()
for v in range(150):
    tmp = b[v].value_counts()
    c = c.append(tmp)
o = c[[b'C', b'G']].sum(axis=1) / c.sum(axis=1) * 100
total = c[[b'C', b'G']].sum(axis=1).sum() / c.sum(axis=1).sum() * 100
o.index = o.index + 1
o.name = 'GC%'
o.to_csv(file2, sep="\t", index_label="pos")
with open(file2, mode="a") as p2:
    p2.write("total\t"+str(total))
