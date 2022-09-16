import numpy as np

with open("data/test.fq") as p:
    d = p.readlines()

i = 0
n = 0
n_sum = 0
for line in d:
    i += 1
    if i % 4 == 0:
        a = np.frombuffer(str(line[0:-1]).encode("utf-8"), dtype=np.uint8)
        n += len(a)
        n_sum += a.sum()

print(n_sum / n)
