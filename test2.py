import numpy as np

with open("data/test.fq") as p:
    d = p.readlines()

i = 0
out = list()
for line in d:
    i += 1
    if i % 4 == 0:
        a = np.frombuffer(str(line[0:-1]).encode("utf-8"), dtype=np.uint8)
        a = (a + 31).tobytes().decode("utf-8")+"\n"
        out.append(a)
        continue
    out.append(line)

with open("data/test.fastq64", mode="w") as p2:
    p2.writelines(out)
