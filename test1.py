
with open("data/test.fq") as p:
    d = p.readlines()

i = 0
out = list()
for line in d:
    i += 1
    if i % 4 == 1:
        tmp1 = ">" + line[1:]
        out.append(tmp1)
    elif i % 4 == 2:
        out.append(line)
with open("data/test.fasta", mode="w") as p2:
    p2.writelines(out)
