#! /usr/bin/python2.6
#! ./
'''
Created on May 7, 2012

@author: Jessica Bonnie
'''


import fix_it
import sys
import intersection
IM_INDEX = 1
BUILD = 'hg18'
MAP_FILE = '/home/jkb4y/cphgdesk_share/Projects/Joe/ImmunoLDBlk16p13.13.map'
LIST_FILE = '/home/jkb4y/cphgdesk_share/Projects/Joe/ListForJessica5-22-2012.txt'
OUT_FILE = '/home/jkb4y/cphgdesk_share/Projects/Joe/KEY_06082012.txt'
INTERSECT_LOC = intersection.LIST_LOC

def old_read_list(list_file):
    line1=True
    im_list = list()
    with open(list_file, mode="r") as listy:
        for line in listy:
            line_list = line.strip().split()
            if line1:
                line1=False
                continue
            elif len(line_list)>1:
                im = line_list[IM_INDEX]
                im_list.append(im)
    return im_list

def read_list(list_loc, annot_dict):
    im_snp_list = list()
    listy = open(list_loc,mode="r")
    for line in listy:
        snp = line.strip()
        im_snp = annot_dict[snp].name
        im_snp_list.append(im_snp)
    return im_snp_list



def look_up_snp(im_name, annot_dict):
    try:
        info =annot_dict[im_name]
        return info
    except KeyError:
        return None




def main():
    build = BUILD
    #list_file = LIST_FILE
    out_file = OUT_FILE
    map_file = MAP_FILE
    intersect_loc = INTERSECT_LOC
    title_list = ['imchip_name','rs_name','hg18_pos','hg19_pos','hg18_lz','intersect']
    
    annot_dict_loc = fix_it.locate_annot_dict(build)
    annot_dict = fix_it.build_annot_dict('MAP',annot_dict_loc)
    intersect_annot_loc = fix_it.locate_annot_dict('hg19')
    intersect_annot_dict = fix_it.build_annot_dict('LOG',intersect_annot_loc)
    
    #im_list = old_read_list(list_file)
    intersect_list = read_list(intersect_loc, intersect_annot_dict)
    index_dict = {'chr':0,'pos':3,'snp':1}
    out = open(out_file, mode="w")
    out.write('\t'.join(title_list)+'\n')
    with open(map_file, mode="r") as mappy:
        for line in mappy:
            line_list = line.strip().split()
            im=line_list[index_dict['snp']]
            info = look_up_snp(im, annot_dict)
            intersect = 'no'

            if info is not None:
                
                if info.name in intersect_list:
                    intersect = 'yes'
                out.write('\t'.join([im,info.rs, info.hg18_pos, info.hg19_pos,
                                     info.lz,intersect]) +'\n')
                
            else:
                out.write('\t'.join([im,im,line_list[index_dict['pos']],
                                     '??','chr16:'+line_list[index_dict['pos']],intersect])+'\n')
    out.close()

if __name__=='__main__':
    main(sys.argv[1:])
