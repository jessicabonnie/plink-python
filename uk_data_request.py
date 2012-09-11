#! /usr/bin/python2.6
#! ./
'''
Created on July 12, 2012

@author: Jessica Bonnie
'''
import sys
import os

def read_list(list_loc):
    snp_list = list()
    with open(list_loc, mode="r") as listy:
        for line in listy:
            rs = line.strip()
            snp_list.append(rs)
    return snp_list

def copy_titles(data_loc, out_loc):
    with open(data_loc, mode = "r") as dat:
        line1 = True
        out = open(out_loc, mode="w")
        for l in dat:
            if line1:
                names = l.strip().split()
                out.write('\t'.join(names)+'\n')
                print l
                line1 = False
    out.close()


def pull_from_data(file_loc, index, snp_list, out_file):
    counter = 0
    list_copy = snp_list[:]
    with open(file_loc, mode="r") as data:
        for line in data:
            lsplit = line.strip().split()
            if lsplit[index] in snp_list:
                out = open(out_file, mode = "a")
                out.write('\t'.join(lsplit)+'\n')
                out.close()
                counter = counter + 1
                list_copy.remove(lsplit[index])
    print("{0} SNPs located in file out of {1}".format(counter, len(snp_list)))
    print list_copy
    

def main(argv):
    base = '/home/jkb4y/work/data/2012Feb1'
    outbase = '/home/jkb4y/work/results/UK_extractions'
    list_loc = os.path.join(outbase,'SNPlist_sent_from_UK.txt')
    cc_loc = os.path.join(base,'hg19','imm_assoc.assoc.logistic')
    cc_out_loc = os.path.join(outbase,'cc_assoc_extract.tbl')
    eurfam_loc = os.path.join(base,'Family_data','eurgdtscan_06062012_lz_hg19.txt')
    eurfam_out_loc = os.path.join(outbase,'eurgdtscan_fam_extract.tbl')
    allfam_loc = os.path.join(base,'Family_data','ALL_fam_plink.tdt')
    allfam_out_loc = os.path.join(outbase,'ALL_fam_extract.tbl')
    eurmeta_loc = os.path.join(base,'eurmeta','eurmeta_06062012_lz_hg19.txt')
    eurmeta_out_loc = os.path.join(outbase,'eurmeta_extract.tbl')
    meta_loc = os.path.join(base,'meta','meta_lz_hg19.tbl')
    meta_out_loc = os.path.join(outbase,'meta_extract.tbl')
    snp_list = read_list(list_loc)
    copy_titles(cc_loc, cc_out_loc)
    copy_titles(eurfam_loc, eurfam_out_loc)
    copy_titles(eurmeta_loc, eurmeta_out_loc)
    copy_titles(meta_loc, meta_out_loc)
    copy_titles(allfam_loc, allfam_out_loc)
    index_dict = dict()
    index_dict['eurfam'] = 0
    index_dict['meta'] = 0
    index_dict['cc']= 1
    index_dict['allfam'] = 1
    pull_from_data(cc_loc, index_dict['cc'], snp_list, cc_out_loc)
    pull_from_data(eurfam_loc,index_dict['eurfam'], snp_list, eurfam_out_loc)
    pull_from_data(allfam_loc,index_dict['allfam'], snp_list, allfam_out_loc)
    pull_from_data(meta_loc, index_dict['meta'], snp_list, meta_out_loc)
    pull_from_data(eurmeta_loc, index_dict['meta'], snp_list, eurmeta_out_loc)

if __name__=='__main__':
    main(sys.argv[1:])
