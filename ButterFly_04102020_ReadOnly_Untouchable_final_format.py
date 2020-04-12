import re
import csv
import operator
from collections import defaultdict
import pandas as pd
import sys
def run():
    count = 0
    totalRow = 0
    # 'Epitope-356214-Folate-Receptor-Alpha.csv'
    if len(sys.argv) != 4:
        print("missing or too many arguments!!!")
        return
    with open(sys.argv[1]) as csvDataFile:
        csvReader = csv.reader(csvDataFile)
        totalRow = sum(1 for _ in csvReader)
        # print('total Row: ' + str(totalRow))

    with open(sys.argv[1]) as csvDataFile:
        csvReader = csv.reader(csvDataFile)
        row1 = next(csvReader)
        count += 1
        file = open('temp.csv', 'w')
        file.write('sample' +',' +'timepoint' + ',' +'percentd_replicate' + ',' +'start'+ ',' + 'end' + ',' + 'Peptide'+ '\n')


        while ('0s' not in row1) and ('.raw' not in row1):
            row1 = next(csvReader)
            count += 1
        cycle = 1
        peptideCheckFlag = int(sys.argv[3])
        while row1 != '':
            # print(row1)
            if cycle == 1:
                peptideCheck = row1[0]
            if peptideCheckFlag == 1:
                #print('*********')
                if peptideCheck != row1[0]:
                    print('++++++++++')
                    peptideCheckFlag = 0
                    print("Wrong Row: " + str(count) + '   ' + row1[2] + '-' + row1[3])
            if 'false' in row1:
                # print('***********************************')
                indexOfBracket = row1[0].find('[')
                peptide = row1[0][:indexOfBracket]
                if indexOfBracket == -1:
                    peptide = row1[0]
                # print(row1[20] + ',' + row1[21] + ',' + row1[35] + ',' + row1[2] + ',' + row1[3] + ',' + peptide )
                file.write(row1[20] + ',' + row1[21] + ',' + row1[35] + ',' + row1[2] + ',' + row1[3] + ',' + peptide +'\n' )
            if count == totalRow:
                break
            row1 = next(csvReader)
            count += 1
            cycle += 1
            if cycle == 1 + int(sys.argv[2]):
                cycle = 1
                    # print(count)

    file.close()
    Ag = ''
    AgAb = ''

    with open('temp.csv') as csvDataFile2:
        csvReader2 = csv.reader(csvDataFile2)
        tempRow = next(csvReader2)
        sortedlist = sorted(csvReader2, key=lambda row: (int(row[3]), int(row[4])), reverse=False)
        dictres = defaultdict(list)
        timepointDic = {}
        peptideList = []
        for x,row in enumerate(sortedlist):
            if x == 0:
                Ag = row[0]
            if x == len(sortedlist) - 1:
                AgAb = row[0]
            timepointDic[row[1]] = timepointDic.get(row[1], 0) + 1
            key = row[0] +'-'+ row[1]
            if row[1] == '0s':
                continue
                # dictres[key] = dictres.get(key, []) + [row[3] + '--' + row[4]]
            elif(row[1] == '15s'):
                peptideList.append(row[5])
                Key0s = row[0] +'-' + '0s'
                dictres[Key0s] = dictres.get(Key0s, []) + [row[3] + '--' + row[4]]
                dictres[key] = dictres.get(key, []) + [row[2]]
            else:
                dictres[key] = dictres.get(key, []) + [row[2]]
        res = dict(dictres)
        keylist  = list(res.keys())
        print('total key: '+ str(len(keylist)))
        print(keylist)

        timepointList = list(timepointDic.keys())
        print(timepointDic)
        header = ''
        header = header + "Peptide"+','+"Peptide Range" + ','
        # for j in range(len(keylist)):
        #     if (j == len(keylist)/2):
        #         continue
        #     if j == 0:
        #         header = header + "Peptide"+','+"Peptide Range" + ','
        #     else:
        #         header = header + keylist[j] +','

        # header = header + ',' + ',' + ',' + 'Peptide Range(Differential)' +','

        
        for k in range(1, len(timepointDic)):
            header = header + timepointList[k] + ','
        
        
        for k in range(1, len(timepointDic)):
            header = header + timepointList[k] + ','

        header += ','


      

        for k in range(1, len(timepointDic)):
            header = header + timepointList[k] + ','


        header = header + ','  + 'Dmax' + ','  +','
        for k in range(1, len(timepointDic)):
            header = header + timepointList[k] + ','

        for k in range(1, len(timepointDic)-1):
            print(timepointDic[timepointList[k]])
            if timepointDic[timepointList[k]] != timepointDic[timepointList[k+1]]:
                return 99

        outputFileName = 'Deuterium Uptake_' + sys.argv[1]
        file2 = open(outputFileName, 'w')
        file2.write(','+','+'D% ' + 'for ' + Ag + ' alone'+','+','+','+','+'D% ' + 'for ' + AgAb+','+','+','+','+','+'Delta D%( ' + Ag+ ' VS '+AgAb +')'+','+','+','+','+','+','+','+'Delta D( ' + Ag+ ' VS '+AgAb +')')


        file2.write('\n')

        file2.write(header)
        file2.write('\n')
        print(keylist)
        for i in range(len(res[keylist[0]])):
            finalres = ''
            finalres += peptideList[2*i]+','
            finalres += res[keylist[0]][i]+','
            for j in range(1,len(keylist)):
                if (j == len(keylist)/2):
                    continue
                if j < len(keylist)//2 :
                    finalres = finalres + "{0:.1f}".format(float(res[keylist[j]][i])) + ','
                else:
                    finalres = finalres  +"{0:.1f}".format(float(res[keylist[j]][i])) + ','
            finalres = finalres +',' 
            DaDiffList = []
            for x in range(1,len(keylist)//2):
                DaDiffList.append(float(res[keylist[x + len(keylist)//2]][i]) - float(res[keylist[x]][i]))
                finalres = finalres + str("{0:.1f}".format(float(res[keylist[x + len(keylist)//2]][i]) - float(res[keylist[x]][i]))) + ','
            finalres += ','
            Dmax = len((peptideList[2*i][2:len(peptideList[2*i])]).replace('P', ''))
            finalres += str(Dmax)
            finalres = finalres + ',' + ',' 
            for p in range(len(DaDiffList)):
                finalres = finalres + str("{0:.1f}".format(DaDiffList[p]*Dmax/100)) +','

            file2.write(finalres)
            file2.write('\n')

  

    file2.close()
    return 0
  



run()