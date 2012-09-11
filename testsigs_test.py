#! /usr/bin/python2.6
#! ./
'''
Created on Jun 28, 2012

@author: Jessica Bonnie
'''
import sys

def check_newthing(single, whole):
    s = open(single, mode="r")
    w = open(whole, mode="r")
    sline1 = True
    wline1 = True
    sdict = dict()
    for sline in s:
        slsplit = sline.strip().split()
        if sline1:
            csnp = slsplit.index('conditional_snp')
            pval = slsplit.index('SNP*_pvalue')
            sline1 = False
        else:
            sdict[slsplit[csnp]] = slsplit[pval]
    s.close()
    for wline in w:
        wlsplit = wline.strip().split()
        if wline1:
            csnp = wlsplit.index('conditional_snp')
            pval = wlsplit.index('SNP*_pvalue')
            wline1 = False
        else:
            if wlsplit[csnp] in sdict:
                if not sdict[wlsplit[csnp]] == wlsplit[pval]:
                    print("OH NO!")
            else:
                print("OH NO!")


def main(argv):
    singlefile = '/home/jkb4y/cphgdesk_share/Projects/IMCHIP/Intersect_SNP_list/TestSigSNPs/Singles/chr17/_17q21.31_z.tbl'
    wholefile = '/home/jkb4y/cphgdesk_share/Projects/IMCHIP/Intersect_SNP_list/TestSigSNPs/chr17/_17q21.31_z.tbl'
    check_newthing(singlefile, wholefile)


if __name__ == '__main__':
    main(sys.argv[1:])
