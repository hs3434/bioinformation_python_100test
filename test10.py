import re
import os

file1 = "data/test.fq"
file2 = "data/test.10id"
key = "CCTGTAATCCCA"
out = list(["", "", "", ""])

if os.path.exists(file2):
    os.remove(file2)

with open(file1) as p, open(file2, mode="a") as p2:
    while 1:
        for v in range(4):
            d = p.readline()
            out[v] = d
        if d == "":
            break
        if re.search(key, out[1]):
            p2.writelines(out[0][:17]+"\n")
