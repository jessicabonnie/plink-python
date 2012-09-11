#! /usr/bin/python2.7
#! ./
'''
Created on Dec 8, 2011

@author: Jessica Bonnie
'''

AA_BIM = '/home/jkb4y/work/Projects/AA/data/AA.bim'
UK_BIM = '/home/jkb4y/work/Projects/UK/data/UK.bim'
RESULT_BIM = '/home/jkb4y/work/Projects/Testing/data/kingtest.bim'

AA_FAM = '/home/jkb4y/work/Projects/AA/data/AA.fam'
UK_FAM = '/home/jkb4y/work/Projects/UK/data/UK.fam'
RESULT_FAM = '/home/jkb4y/work/Projects/Testing/data/kingtest.fam'

AA_BED = '/home/jkb4y/work/Projects/AA/data/AA.bed'
UK_BED = '/home/jkb4y/work/Projects/UK/data/UK.bed'
RESULT_BED = '/home/jkb4y/work/Projects/Testing/data/kingtest.bed'

def make_test_bim(AA_bim, UK_bim, result):
    counter = 1
    out = open(result, mode="w")
    with open(AA_bim, mode="r") as AA:
        for line in AA:
            if counter%2==0:
                line_list = line.split()
                line_list[3] = str(int(line_list[3]) + 1)
                out.write(line + '\t'.join(line_list) + '\n')
            else:
                out.write(line)
            counter +=1
    counter = 1
    with open(UK_bim, mode = "r") as UK:
        UK.next()
        for line in UK:
            if counter < 1000:
                out.write(line + '\n')
                counter +=1
    out.close()

def make_test_fam(AA_fam, UK_fam, result):
    counter = 1
    out = open(result, mode="w")
    with open(AA_fam, mode="r") as AA:
        for line in AA:
            if counter%2==0:
                line_list = line.split()
                line_list[1] = str(int(line_list[1]) + 1)
                out.write(line + ' '.join(line_list) + '\n')
            else:
                out.write(line)
            counter +=1
    counter = 1
    with open(UK_fam, mode = "r") as UK:
        UK.next()
        for line in UK:
            if counter < 1000:
                out.write(line)
                counter +=1
    out.close()

def make_test_bed(AA_bed, UK_bed, result):
    counter = 1
    out = open(result, mode="w")
    with open(AA_bed, mode="r") as AA:
        for line in AA:
            if counter%2==0:
                out.write(line + line)
            else:
                out.write(line)
            counter +=1
    counter = 1
    with open(UK_bed, mode = "r") as UK:
        for line in UK:
            if counter < 1000:
                out.write(line)
                counter +=1
    out.close()
                
def main():
    out = open(RESULT_BIM, mode = "w")
    with open(AA_BIM, mode = "r") as AA:
        for line in AA:
            out.write(line)
    out.close()
    #make_test_bim(AA_BIM,UK_BIM, RESULT_BIM)
    make_test_fam(AA_FAM,UK_FAM, RESULT_FAM)
    make_test_bed(AA_BED, UK_BED, RESULT_BED)

if __name__=='__main__':
    main()
