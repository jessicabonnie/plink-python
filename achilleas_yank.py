#! /usr/bin/python2.6
#! ./
'''
Created Oct 28, 2011

@author: Jessica Bonnie
'''


import os
import re
import sys
import getopt
import subprocess

import pc_toolbox
import fix_it
import meta_yank

BUILD = 'hg19'
REGION_LOC = '/home/jkb4y/work/data/Region_Lists/hg19/T1D_regions_hg19_dbase_and_sig.txt'
ACHILLEAS_KEY_LOC ='/home/jkb4y/cphgdesk_share/Achilleas/cis-eQTLs/immIDs2currentIDs.tsv'
AP_SNP_TITLE = 'SNP_AP'
IM_SNP_TITLE = 'SNP_IM'
JB_SNP_TITLE = 'SNP_LZ_JB'
LZ_TITLE ='LZ' #this is the lz col in case a kicked snp needs to be drawn
EXTREME_LZ = 'EXTREME_LZ'#lz "make it work" col, all snps in chr:pos form to make sure lz can "see" them
ANNOT_TITLE = 'annotation'



P_TITLE = 'P'

def create_eqtl_index(table_loc, perm=False):
    line1=True
    with open(table_loc, mode="r") as table:
        for line in table:
            if line1:
                if perm:
                    index_dict = pc_toolbox.read_perm_titles(line)
                else:
                    index_dict = pc_toolbox.read_eqtl_titles(line)
                return index_dict

def read_eqtl_titles(title_line, perm=False):
    title_list = title_line.strip().split()
    a1_index = None
    ci_hi_index = None
    ci_lo_index = None
    or_index = None
    stat_index = None
    lz_choice = JB_SNP_TITLE
    print lz_choice
    print JB_SNP_TITLE
    if perm:
        p_name = 'EMP1'
        p_index =  title_list.index(p_name)
    if not perm:
        p_name = 'P'
        a1_index = title_list.index('A1')
        or_index = title_list.index('BETA')
        stat_index = title_list.index('STAT')
        p_index = title_list.index(p_name)
    else:
        p_index =  title_list.index('EMP1')
        
    index_dict = {'p':p_index,
                  'chr':title_list.index('CHR'),
                  'pos':title_list.index('BP'),
                  'snp':title_list.index(lz_choice),
                  'im':title_list.index(IM_SNP_TITLE),
                  'or':or_index,
                  't':stat_index,
                  'hi':ci_hi_index,
                  'lo':ci_lo_index,
                  'a1':a1_index,
                  'pvalcol':p_name,
                  'markercol':lz_choice}
    return index_dict

def build_key(key_loc, purpose = 'IM'):
    key_dict = dict()
    with open(key_loc, mode="r")as key:
        for line in key:
            lsplit = line.strip().split()
            if purpose == 'IM':
                key_dict[lsplit[1]]=lsplit[0]
            if purpose == 'AP':
                key_dict[lsplit[0]]=lsplit[1]
    return key_dict

def write_multi_table(table_loc, multi_loc, key_dict, annot_dict):
    new = open(multi_loc, mode="w")
    base, ext = os.path.splitext(table_loc)
    kicked_loc = base + '_kick'+ext
    kicked = open(kicked_loc, mode="w")
    snp_index = 2
    chr_index = 3
    bp_index = 4
    line1 = True
    with open(table_loc, mode="r") as old:
        for line in old:
            lsplit = line.strip().split()
            if line1:
                snp_index = lsplit.index('SNP')
                chr_index = lsplit.index('CHR')
                bp_index = lsplit.index('BP')
                lsplit[snp_index]=AP_SNP_TITLE
                lsplit.append(IM_SNP_TITLE)
                kicked.write('\t'.join(lsplit)+'\n')
                lsplit.extend([JB_SNP_TITLE,LZ_TITLE,ANNOT_TITLE,EXTREME_LZ])
                new.write('\t'.join(lsplit)+'\n')
                line1 = False
            else:
                snp = lsplit[snp_index]
                im = key_dict[snp]
                lsplit.append(im)
##                chrom = lsplit[chr_index]
##                pos = lsplit[bp_index]
##                elz = "chr{0}:{1}".format(chrom,bp)
##                lsplit.append(elz)
                
                try:
                    annot = fix_it.get_lz_annot(annot_dict,im)
                    lz = annot_dict[im].lz
                    chrom = annot_dict[im].hg19_chr
                    pos = annot_dict[im].hg19_pos
                    elz = "chr{0}:{1}".format(chrom,pos)
                    lsplit.extend([lz,lz,annot,elz])
##                    lsplit.append(annot)
##                    lsplit.append(elz)
                except KeyError:
                    kicked.write('\t'.join(lsplit)+'\n')
                    lz_jb = 'NA'
                    lz = snp
                    chrom = lsplit[chr_index]
                    pos = lsplit[bp_index]
                    elz = "chr{0}:{1}".format(chrom,pos)
                    annot = 'no-annotation'
                    lsplit.extend([lz_jb,lz,annot,elz])
                    
                new.write('\t'.join(lsplit)+'\n')
    kicked.close()
    new.close()
            

def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global table_loc, region_loc, out_base, annot, build, fix_args
    global weight_min, freq_loc, covar_loc, bfile
    global family
    
    table_loc = None
    region_loc = REGION_LOC
    out_base = None
    freq_loc = None
    weight_min = None#WEIGHT_MIN
    sex_weight = None#SEX_WEIGHT
    covar_loc = None
    bfile = None
    annot=None
    build = BUILD
    fix_args = list()
    family = False
    
    try: 
        opts, args = getopt.getopt(argv, "ht:r:o:f:b:",
                                   ["help","table=","region-list=",
                                    "out=","freq=","bfile=",
                                    "covar=","annot=","build=","family="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-f","--freq"):
            freq_loc = arg
        elif opt in ("-t","--table"):
            table_loc = arg
        elif opt in ("-r","--region-list"):
            region_loc = arg
        elif opt in ("-o","--out"):
            out_base = arg
        elif opt in ("-b","--bfile"):
            bfile = arg
        elif opt in ("--family"):
            meta_loc = arg
            family = True
        elif opt in ("--covar"):
            covar_loc = arg
        elif opt in ("--build"):
            build = arg
            fix_args.extend([opt,arg])
##        elif opt in ("--fix"):
##            fix = True
        elif opt in ("--annot"):
            annot = arg




def main(argv):
    global table_loc, region_loc, achilleas_key, build,out_base
    weight_min = 0
    cl_arguments(argv)
    perm = False
    if 'perm' in table_loc:
        perm = True
    achilleas_key_loc = ACHILLEAS_KEY_LOC
    achilleas_dict = build_key(achilleas_key_loc, 'IM')
    annot_dict_loc = fix_it.locate_annot_dict(build)
    annot_dict = fix_it.build_annot_dict('MAP',annot_dict_loc)
    base, ext = os.path.splitext(table_loc)
    multi_loc = base + '_JB'+ext
    write_multi_table(table_loc, multi_loc, achilleas_dict,annot_dict)
    
##    region_list = pc_toolbox.create_region_list(region_loc)
##    keep_loc = out_base+'_plink_keep.txt'
##    plink_out = out_base+'_sigs_eQTL'
    index_dict = create_eqtl_index(multi_loc, perm)
    log_dict = fix_it.build_annot_dict('LOG',annot_dict_loc)
##    z_list, lw_list = meta_yank.find_leastP_SNPs(region_list, multi_loc,
##                                       index_dict, out_base, freq_loc,
##                                       keep_loc, weight_min, 'eqtl',log_dict)
##    print z_list




if __name__ == '__main__':
    main(sys.argv[1:])




