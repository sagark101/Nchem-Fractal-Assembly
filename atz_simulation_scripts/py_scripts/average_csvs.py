import sys

with open(sys.argv[1], 'r') as csvfiles:
    list_of_csvs = csvfiles.readlines()

final = []
for csvf in list_of_csvs:
    csvf = csvf.rstrip('\r\n')
    with open(csvf, 'r') as csv:
        mylist = csv.readlines()
    
    count = 0
    master = mylist[0].split(',')
    master[1] = 0
    master[-1] = 0#float(master[-1])
    master[-3] = 0#float(master[-3])
    master[-4] = 0#float(master[-4])
    master[-2] = 0
    for line in mylist:
        
        listL = line.split(',')
        if int(listL[-3]) + int(listL[-4]) >= 5000:
            master[1] += float(listL[1])
            master[-1] += float(listL[-1])
            master[-3] += float(listL[-3])
            master[-4] += float(listL[-4])
            if float(listL[-2]) > 3.0:
                master[-2] += 3.0
            else:
                master[-2] += float(listL[-2])
            count += 1
    
    if count > 0:
        master[1] = master[1]/float(count)
        master[-1] = master[-1]/float(count)
        master[-2] = master[-2]/float(count)
        master[-3] = master[-3]/float(count)
        master[-4] = master[-4]/float(count)

    master = [str(x) for x in master]
    master_str = str(count) + ',' + ','.join(master) + '\n'
    final.append(master_str)

with open('master_frac_null_temp9.csv', 'w') as mycsv:
    mycsv.writelines(final)
