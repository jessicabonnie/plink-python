#! /ust/bin/python2.7
#! ./
'''
Created Jan 25, 2012

@author: Jessica Bonnie
'''
from collections import namedtuple
import os
import sys
import getopt
import shutil
import fix_it_fam
import fix_it
import meta_yank
import pc_toolbox


import table_yank


FAM_CHR_TITLE = 'CHR'
FAM_SNP_TITLE = 'SNP'
FAM_P_TITLE = 'PVALUE'
FAM_POSITION_TITLE = 'POS'
FAM_Z_TITLE = 'Z_GDT'
FAM_WEIGHT_TITLE = 'Pair'
FAM_A1_TITLE = 'AL1'
FAM_A2_TITLE = 'AL2'



def read_fam_titles(line):
    title_list = line.strip().split()
    meta_chr_col = title_list.index(FAM_CHR_TITLE)
    p_col = title_list.index(FAM_P_TITLE)
    snp_col = title_list.index(FAM_SNP_TITLE)
    position_col = title_list.index(FAM_POSITION_TITLE)
    z_col = title_list.index(FAM_Z_TITLE)
    weight = title_list.index(FAM_WEIGHT_TITLE)
    meta_a = title_list.index(FAM_A1_TITLE)
    a2 = title_list.index(FAM_A2_TITLE)
    fam_indices = {'chr':meta_chr_col, 'p':p_col,
                    'snp':snp_col,'z':z_col,
                    'pos':position_col,'weight':weight,
                    'a1':meta_a, 'a2':a2}
    #print("Index Dict:")
    #print fam_indices
    return fam_indices


def retrieve_fam_data(snp_list, fam_loc, out_loc, annot_dict):
    line1 = True
    out = open(out_loc, mode = "w")
    with open(fam_loc, mode = "r") as fam:
        for line in fam:
            if line1:
                index_dict = pc_toolbox.read_fam_titles(line)
                title_list = ['lz_snp','rs_snp','imchip_snp',FAM_CHR_TITLE,FAM_POSITION_TITLE,
                              FAM_P_TITLE, FAM_Z_TITLE,FAM_WEIGHT_TITLE,FAM_A1_TITLE,FAM_A2_TITLE]
                out.write('\t'.join(title_list)+'\n')
                line1 = False
            else:
                line_list = line.strip().split()
                snp = line_list[index_dict['snp']]
                weight = line_list[index_dict['weight']]
                z = line_list[index_dict['z']]
                pos = line_list[index_dict['pos']]
                a1 = line_list[index_dict['a1']]
                a2 = line_list[index_dict['a2']]
                chro = line_list[index_dict['chr']]
                pval = line_list[index_dict['p']]
                if snp in snp_list:
                    rs = annot_dict[snp].rs
                    im = annot_dict[snp].name
                    new_list = [snp,rs, im,chro,pos,pval,z,weight,a1,a2]
                    #or annot_dict[snp].name in snp_list or annot_dict[snp].rs in snp_list:
                    out.write('\t'.join(new_list)+'\n')
    out.close()

def get_snp_list(meta_yank_loc):
    line1 = True
    snp_list = list()
    counter = 0
    with open(meta_yank_loc, mode = "r") as meta:
        for line in meta:
            if line1:
                yank_indices = pc_toolbox.read_yank_titles(line)
                line1 = False
            else:
                line_list = line.strip().split()
                snp_list.append(line_list[yank_indices['lz']])
                counter = counter + 1
    print("I find {0} SNPs in the meta_yank list.".format(counter))
    return snp_list
                
            

def main(argv):
    meta_yank_loc = '/home/jkb4y/work/results/Intersection/eurmeta_06062012/RegionYank/eurmeta_yank.tbl'
    fam_loc = '/home/jkb4y/work/data/2012Feb1/Family_data/eurgdtscan_06062012_lz_hg19_intersect.txt'
    out_loc = '/home/jkb4y/work/results/Intersection/family/eurmeta_snp_info.tbl'
    yank_out = '/home/jkb4y/work/results/Intersection/family/RegionYank/family'
    BFILE = '/home/jkb4y/work/data/2012Feb1/intersect_fam_uk/hg19/intersect'
    FREQ = '/home/jkb4y/work/data/2012Feb1/intersect_fam_uk/hg19/intersect_controls.frq'
    REGION_LOC = '/home/jkb4y/work/data/Region_Lists/hg19/T1D_regions_hg19_05242012.txt'
    COVAR = '/home/jkb4y/work/data/2012Feb1/intersect_fam_uk/hg19/intersect.cov'
    build = 'hg19'
##
##    table_yank.main(['--family',fam_loc, '--build',build,
##                     '--bfile',BFILE,'--freq',FREQ,'-r',REGION_LOC,
##                     '--covar',COVAR, '-o', yank_out])
    annot_dict_loc = fix_it.locate_annot_dict(build)
    annot_dict = fix_it.build_annot_dict('LOG',annot_dict_loc)
    snp_list = get_snp_list(meta_yank_loc)
    fixed_table_loc = fam_loc
    #fixed_table_loc = fix_it.locate_fixed_table(fam_loc, build)
    retrieve_fam_data(snp_list, fixed_table_loc, out_loc, annot_dict)
##
    meta_yank.main(['--family',fam_loc, '--build',build,
                     '--bfile',BFILE,'--freq',FREQ,'-r',REGION_LOC,
                     '--covar',COVAR, '-o', yank_out])
    

if __name__ == '__main__':
    main(sys.argv[1:])
    
