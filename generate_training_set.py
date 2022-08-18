"""
This script takes the scores.dat file generated by MCC.py and generates a table with the following information:
| name | sequence | correct structure | ExpertRNA structure | MCC score of predicted vs correct |
"""

from sys import argv
from json import dump

output_data = {}

# argv 1 is the scores.dat file
with open(argv[1]) as f:
    data = f.readlines()
    for line in data:
        l = line.split(' ')
        output_data[l[0]] = {'MCC_expert' : float(l[2])}
        output_data[l[0]]['MCC_RNAfold'] = float(l[1])

#argv 2 is the Results file
with open(argv[2]) as f:
    data = f.readlines()[1::4]
    for line in data:
        l = line.split(',')
        n = l[1]
        if n in output_data:
            output_data[n]['seq'] = l[2]
            output_data[n]['corr_str'] = l[3]
            output_data[n]['RNAfold_str'] = l[4]
            output_data[n]['pred_str'] = l[5]

#argv 3 is the output file
with open(argv[3], 'w+') as f:
    #dump(output_data, f)
    f.write("name,sequence,correct_structure,RNAfold_structure,expertRNA_structure,MCC_RNAfold,MCC_expert\n")
    for k in output_data.keys():
        f.write("{},{},{},{},{},{},{}\n".format(k, output_data[k]['seq'], output_data[k]['corr_str'], output_data[k]['RNAfold_str'], output_data[k]['pred_str'], output_data[k]['MCC_RNAfold'], output_data[k]['MCC_expert']))