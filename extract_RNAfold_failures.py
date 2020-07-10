# Run after names are fixed and failures are extracted into a file.
#failures are extracted with the following bash command
#cat scores_1.dat | awk '{if ($2 < 0.1) {print}}' | cut -d ' ' -f 1 > test.txt

from os import listdir

fnames = []
with open('test.txt', 'r') as f: 
    data = f.readlines()
    for d in data:
        fnames.append(d.strip())

out = ''
files = listdir()
for fi in files:
    if 'fixed' in fi:
        with open(fi, 'r') as f:
            data = f.readlines()
            for d in data:
                if d.split(',')[0] in fnames:
                    fnames.remove(d.split(',')[0])
                    out = out + d

with open('bad_RNAfold.csv', 'w+') as f:
    f.write(out)