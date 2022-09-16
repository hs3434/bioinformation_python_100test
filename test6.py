import numpy as np

file1 = "data/test.fq"
file2 = "data/test.total"
a = np.empty((0, 150))
n = 0
with open(file1) as p:
    while 1:
        d = p.readline()
        if d == "":
            break
        n += 1
        if n % 4 == 0:
            tmp = np.frombuffer(str(d[0:-1]).encode("utf-8"), dtype=np.uint8)
            a = np.row_stack((a, tmp))
a2 = a.mean(0)
i = a.shape[0] // 2
a.sort(0)
with open(file2, mode="w") as p2:
    p2.write("pos\tmean\tmedian\n")
    for v in range(150):
        p2.write(str(v + 1) + "\t" + str(a2[v]) + "\t" + str(a[i, v])+"\n")
