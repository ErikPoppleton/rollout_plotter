#Run this on new datasets coming in from Rollout.  
#There are frequently duplicate names with different sequences
#Need to split out the names

from os import listdir
files = listdir()

names = {}
previous_name = ''
sub_name_count = 0
out = ''
for fi in files:
    if 'Result' in fi:
        with open(fi, 'r') as f:
            data = f.readlines()
            for d in data:
                name = d.split(',')[0]
                if name != previous_name or sub_name_count > 4:
                    sub_name_count = 0

                if name in names.keys() and name != 'NONE' and sub_name_count == 0:
                    names[name] += 1
                    d = d.replace(name, name+'_'+str(names[name]))
                    sub_name_count += 1
                elif name in names.keys() and name != 'NONE':
                    d = d.replace(name, name+'_'+str(names[name]))
                    sub_name_count += 1
                elif name != 'NONE':
                    names[name] = 0
                    d = d.replace(name, name+'_'+str(names[name]))
                    sub_name_count += 1

                previous_name = name
                out += d

with open ('Res_fixed.csv', 'w') as f:
    f.write(out)
