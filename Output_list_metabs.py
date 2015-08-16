#Need the counts dict of metab to generated the list of metabs.
li = counts.items()
li.sort()
li.sort(key=lambda x: x[1])
a = []
b = []
for i in li:
    #a.append(list(i))
    b.append(str(i[0]) + "/" + str(i[1]))
with open("conteo_metabs_.txt", "w") as lista:
    for x in b:
        lista.write(x+"\n")
    lista.close()
#Could be used in recon and HMRA files.
