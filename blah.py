#! /usr/bin/python2.7

import dup_check
import duplicate_check

def main():
    chromsome = 2
    chrband = '2q33.2'
    #merge_out = '/home/jkb4y/work/results/HapMap/NIH_HapMap/chr{0}/{1}'.format(chromsome,chrband)
    merge_out = '/home/jkb4y/work/results/HapMap/NIH_HapMap/PLINK/all'.format(chromsome,chrband)
    dup_check.summarize_diff(merge_out)
    duplicate_check.summarize_diff(merge_out)

if __name__=='__main__':
    main()
    
