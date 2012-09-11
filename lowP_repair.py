#! /ust/bin/python2.6
#! ./
'''
Created Mar 6, 2012

@author: Jessica Bonnie
'''

from collections import namedtuple
import os
import sys


import meta_toolbox
import pc_toolbox
import fix_it

EURMETA = '/home/jkb4y/work/data/2012Feb1/eurmeta/eurmeta.tbl'
EURLOW = '/home/jkb4y/work/data/2012Feb1/eurmeta/eurlowPdata_JBedit.txt'
OUT_LOC = '/home/jkb4y/work/data/2012Feb1/eurmeta/eurmeta_06062012.txt'

def read_lowP(lowP_table):
    line1 = True
    fix_dict = dict()
    LineFo = namedtuple("LineFo","snp,a1,a2,w,z,p,direction,chro,pos")
    with open(lowP_table, mode="r") as lowP:
        for line in lowP:
            if line1:
                line1 = False
                continue
            else:
                line_split = line.strip().split()
                chro = line_split[6]
                pos = line_split[3]
                a1 = line_split[10]
                a2 = line_split[11]
                snp = line_split[7]
                z = line_split[17]
                weight = line_split[18]
                direction = line_split[21]
                line_fo = LineFo(snp=snp, a1=a1,a2=a2,
                                 w=weight,z=z,p='0',
                                 direction=direction,
                                 chro=chro,pos=pos)
                fix_dict[snp]= line_fo
    return fix_dict
def read_lowPJB(lowP_table):
    line1 = True
    fix_dict = dict()
    LineFo = namedtuple("LineFo","snp,a1,w,z")
    with open(lowP_table, mode="r") as lowP:
        for line in lowP:
            if line1:
                line1 = False
                continue
            else:
                line_split = line.strip().split()
                a1 = line_split[1]
                snp = line_split[0]
                z = line_split[7]
                weight = line_split[8]
                line_fo = LineFo(snp=snp, a1=a1,
                                 w=weight,z=z)
                fix_dict[snp]= line_fo
    return fix_dict


def fix_eur(orig_table, out_table, fix_dict):
    line1 = True
    out = open(out_table, mode="w")
    key_list = fix_dict.keys()
    with open(orig_table, mode="r") as orig:
        for line in orig:
            if line1:
                index_dict = pc_toolbox.read_meta_titles(line)
                ls = line.strip().split()
                ls.append('eurlowP_SNP')
                out.write('\t'.join(ls)+'\n')
                #out.write(line)
                line1 = False
            else:
                lp_flag = '.'
                lsplit = line.strip().split()
                snp = lsplit[index_dict['snp']]
                if snp in key_list:
                    lsplit[index_dict['z']]=fix_dict[snp].z
                    #lsplit[index_dict['p']]=fix_dict[snp].p
                    lsplit[index_dict['weight']]=fix_dict[snp].w
                    lsplit[index_dict['a1']]=fix_dict[snp].a1
                    lsplit[index_dict['a2']]='.'
                    lsplit[index_dict['direction']]='.'
                    key_list.remove(snp)
                    lp_flag = 'yes'
                lsplit.append(lp_flag)
                out.write('\t'.join(lsplit) + '\n')
    if len(key_list) > 0:
        print("{0} keys not found.".format(len(key_list)))
        annot_dict_loc = fix_it.locate_annot_dict('hg18')
        annot_dict = fix_it.build_annot_dict('MAP', annot_dict_loc)
        for key in key_list:
            #print fix_dict
            print(key)
            a1 = fix_dict[key].a1
            a2 = '.'
            weight = fix_dict[key].w
            z = fix_dict[key].z
            p = '0'
            direction = '.'
            chro = annot_dict[key].hg18_chr
            pos = annot_dict[key].hg18_pos
            lp_flag = 'yes'
            key_info = [key,a1,a2,weight,z,p,direction,chro,pos,lp_flag]
            out.write('\t'.join(key_info)+'\n')
            #out.write('\t'.join(fix_dict[key])+'\n')
    out.close()

def main(argv):
                         
    fix_dict = read_lowPJB(EURLOW)
    fix_eur(EURMETA, OUT_LOC, fix_dict)
    
if __name__=='__main__':
    main(sys.argv[1:])
