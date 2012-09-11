#! /usr/bin/python2.7
'''
Created on Sep 14, 2011

@author: Jessica Bonnie
'''
import sys
import duplicate_check

NOVEMBER = '/home/jkb4y/work/data/2011Nov2/copy_of_UK.bim'
HG18= '/home/jkb4y/work/data/2012Feb1/hg18/imm.bim'
HG19= '/home/jkb4y/work/data/2012Feb1/hg19/hg19.bim'
OUT = '/home/jkb4y/work/data/2012Feb1/CHECKME.txt'
SNP_SAME = '/home/jkb4y/work/data/2012Feb1/NAME_SAME_POS_DIFF.txt'
POS_SAME = '/home/jkb4y/work/data/2012Feb1/POS_SAME_NAME_DIFF.txt'
ONLY_ONE = '/home/jkb4y/work/data/2012Feb1/ONLY_IN_ONE.txt'
def check_positions():
    hg18_list = duplicate_check.list_info(HG18)
    nov_list = duplicate_check.list_info(NOVEMBER)
    results = open(OUT, mode = "w")
    snp_same = open(SNP_SAME, mode = "w")
    pos_same = open(POS_SAME, mode = "w")
    only_one = open(ONLY_ONE, mode = "w")
    title_list = ['nov_snp','nov_chr','nov_pos',
                  'hg18_snp','hg18_chr','hg18_pos']
    results.write('\t'.join(title_list)+'\n')
    snp_same.write('\t'.join(title_list)+'\n')
    pos_same.write('\t'.join(title_list)+'\n')
    only_one.write('\t'.join(title_list)+'\n')
    blank = ['--','--','--']
    listy=list()
    listy2 = list()
    for fo in hg18_list:
        if fo not in nov_list:
            for nov in nov_list:
                if fo.snp == nov.snp:
                    print("SNP names same!")
                    print nov
                    print fo
                    results.write('\t'.join(nov)+'\t'+'\t'.join(fo)+'\n')
                    snp_same.write('\t'.join(nov)+'\t'+'\t'.join(fo)+'\n')
                    listy.append(fo)
                    listy.append(nov)
                if nov.chro == fo.chro and nov.pos == fo.pos:
                    print("SNP positions same!")
                    print nov
                    print fo
                    results.write('\t'.join(nov)+'\t'+'\t'.join(fo)+'\n')
                    pos_same.write('\t'.join(nov)+'\t'+'\t'.join(fo)+'\n')
                    listy.append(fo)
                    listy.append(nov)
            if fo not in listy:
                results.write('\t'.join(blank)+'\t'+'\t'.join(fo)+'\n')
                only_one.write('\t'.join(blank)+'\t'+'\t'.join(fo)+'\n')
##    listy2=list()
    for nov in nov_list:
        if nov not in hg18_list:
            for fo in hg18_list:
                if nov.snp == fo.snp:
                    print("SNP names same!")
                    print nov
                    print fo
                    results.write('\t'.join(nov)+'\t'+'\t'.join(fo)+'\n')
                    snp_same.write('\t'.join(nov)+'\t'+'\t'.join(fo)+'\n')
                    listy.append(nov.snp)
                    listy.append(fo.snp)
                if nov.chro == fo.chro and nov.pos == fo.pos:
                    print("SNP positions same!")
                    print nov
                    print fo
                    results.write('\t'.join(nov)+'\t'+'\t'.join(fo)+'\n')
                    pos_same.write('\t'.join(nov)+'\t'+'\t'.join(fo)+'\n')
                    listy.append(nov.snp)
                    listy.append(fo.snp)
            if nov not in listy:
                results.write('\t'.join(nov)+'\t'+'\t'.join(blank)+'\n')
                only_one.write('\t'.join(nov)+'\t'+'\t'.join(blank)+'\n')
    print listy
    print listy
    results.close()
    only_one.close()
    pos_same.close()
    snp_same.close()
            



def main(argv):
    check_positions()

if __name__=='__main__':
    main(sys.argv[1:])
