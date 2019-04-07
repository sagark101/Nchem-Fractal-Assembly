import sys

with open(sys.argv[1], 'r') as myFile:
    summaryfile = myfile.readlines()

final_dictionary = {}

for null in [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
    for frac in [1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]:
        key = 'null'+str(null) + '_' + 'frac'+str(frac)
        #dimension, lacunarity, member count
        final_dictionary[key] = [ 0.0, 0.0, 0.0 ]

for data_line in summaryfile:
    split_line = data_line.split(',')
    keytype = '_'.join(split_line[0].split('_')[4:6])
    final_dictionary[keytype][0] += float(split_line[1])
    final_dictionary[keytype][1] += float(split_line[2])
    final_dictionary[keytype][2] += 1

outlines = []
for k,v in final_dictionary.iteritems():
    averageD = v[0]/v[2]
    averageL = v[1]/v[2]
    outlines.append(k + ',' + str(averageD) + ',' + str(averageL) + '\n')

with open(sys.argv[2], 'w') as newf:
    newf.writelines(outlines)

