#! /usr/bin/python2.6
#! ./
'''
Created Oct 28, 2011

@author: Jessica Bonnie
'''
import sys
import os
import pc_toolbox
import fix_it

CC_LOC = '/home/jkb4y/work/data/2012Feb1/hg19/imm.bim'
FAM_LOC= '/home/jkb4y/work/data/2012Feb1/Family_data/hg19/eur.bim'
FAM_TABLE = '/home/jkb4y/work/data/2012Feb1/Family_data/eurgdtscan_08032012_lz_hg19.txt'
LIST_LOC = '/home/jkb4y/work/data/2012Feb1/intersect_hg19.list'
META_LOC = '/home/jkb4y/work/data/2012Feb1/eurmeta/eurmeta_08032012_lz_hg19.txt'
CHR_TABLE = '/home/jkb4y/work/data/2012Feb1/eurmeta/chromsome_table_hg19_16.txt'
BASE_LIST_LOC = '/home/jkb4y/work/data/2012Feb1/intersect.list'
def read_map(map_loc):
    listy = list()
    with open(map_loc, mode="r")as mappy:
        for line in mappy:
            line_split = line.strip().split()
            listy.append(line_split[1])
    return listy

def make_list(cc_loc, fam_loc, list_loc, annot_dict):

    cc_list = read_map(cc_loc)
    fam_list = read_map(fam_loc)
    intersect_list = list()
    base, ext = os.path.splitext(list_loc)
    hg19_list = base + '_hg19.list'
    im_list = base + '.list'
    im_out = open(im_list, mode = "w")
    out = open(hg19_list,mode="w")
    for snp in cc_list:
        if snp in fam_list:
            intersect_list.append(snp)
            out.write(snp +'\n')
            im = annot_dict[snp].name
            im_out.write(im + '\n')
    out.close()
    im_out.close()
    return intersect_list

def read_list(intersect_loc):
    intersect_list = list()
    with open(intersect_loc, mode="r") as listy:
        for line in listy:
            strip = line.strip()
            intersect_list.append(strip)
    return intersect_list

def filter_table(meta_loc, list_loc,table_type='META'):
    intersect_list = read_list(list_loc)
    base, ext = os.path.splitext(meta_loc)
    new_meta_loc = base + '_intersect'+ext
    line1 = True
    meta = open(meta_loc, mode="r")
    new_meta = open(new_meta_loc, mode="w")
    for line in meta:
        if line1:
            if table_type == 'META':
                index_dict = pc_toolbox.read_meta_titles(line)
            elif table_type == 'FAMILY':
                index_dict = pc_toolbox.read_fam_titles(line)
            new_meta.write(line)
            line1 = False
        else:
            line_list = line.strip().split()
            snp = line_list[index_dict['snp']]
            for isnp in intersect_list:
                if snp == isnp:
                    #print("Found in both lists: {0}".format(snp))
                    new_meta.write(line)
                    continue
    meta.close()
    new_meta.close()
                        

def main(argv):
    cc_loc = CC_LOC
    fam_loc = FAM_LOC
    list_loc = LIST_LOC
    meta_loc = META_LOC
    chr_table = CHR_TABLE
    annot_loc = fix_it.locate_annot_dict('hg19')
    annot_dict = fix_it.build_annot_dict('LOG',annot_loc)
    intersect_list = make_list(cc_loc, fam_loc, BASE_LIST_LOC, annot_dict)
    #intersect_list = read_list(list_loc)
    filter_table(meta_loc, list_loc,'META')
    #intersect_list = read_list(list_loc)
    #filter_table(chr_table, list_loc)
    filter_table(FAM_TABLE, list_loc,'FAMILY')

    

if __name__ == '__main__':
    main(sys.argv[1:])
    
