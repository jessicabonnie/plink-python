#! /usr/bin/python2.7
'''
Created on Sep 14, 2011

@author: Jessica Bonnie
'''

import os
import sys
import getopt
import shutil
from collections import namedtuple

import pc_toolbox

MAP_LOC = None
DUP_LOC = None
REPAIR_LOC = None

#REPAIR_DICT = {'imm_2_230785920':'rs11556887','imm_2_230857352':'rs4972946'}

def fix_map(map_loc, repair_loc):
    '''
    Replaces all SNP names that do not begin with 'rs' with a name using
    chr<chromosome#>:<basepair position> format.
    Args:
        map_loc -- filepath of bim or map file
    
    '''
    (basepath, ext) = os.path.splitext(map_loc)
    orig_rename = basepath + '~'
    new_loc = str(rename_as_necessary(orig_rename, ext)) + ext
    shutil.copy(map_loc, new_loc)
    repair_dict = pc_toolbox.read_dict(repair_loc)
    
    with open(new_loc, mode='r') as original_map:
        with open (map_loc, mode='w') as new_map:
            for orig_line in original_map:
                orig_list = orig_line.split()
##                if orig_list[1] in dup_list:
##                    if orig_list[1] in repair_dict:
##                        new_snp = repair_dict[orig_list[1]]
##                        orig_list[1] = new_snp
##                else:
##                    new_snp = 'chr'+ orig_list[0]+':'+ orig_list[3]
##                    orig_list[1] = new_snp
                if orig_list[1].startswith('rs'):
                    new_snp = orig_list[1]
                elif orig_list[1] in repair_dict:
                    new_snp = repair_dict[orig_list[1]]
                else:
                    new_snp = 'chr'+ orig_list[0]+':'+ orig_list[3]
##                    if orig_list[1] in repair_dict:
##                        new_snp = repair_dict[orig_list[1]]
##                    else:
##                        new_snp = 'chr'+ orig_list[0]+':'+ orig_list[3]
                orig_list[1] = new_snp
                new_map.write('\t'.join(orig_list)+'\n')                

                
def rename_as_necessary(new_orig,ext):
    '''
    Determines if the desired new name for the input map file already exists,
    and, if so, chooses another name.
    Args:
        new_orig -- target to which to map file could be moved
        ext -- extension of the map file
    Returns:
        next_rename -- target name to which the original map file can
                        be saved without overwriting any other file
    '''
    next_rename = new_orig
    if os.path.exists(new_orig + ext):
        print( '''There is already a file named {0}, which means that map_adapt has
already been run on this map file. To avoid overwriting the contents of {0},
the contents of the input file will be written to: '''.format(new_orig + ext))        
        if new_orig[-1].isdigit():
            pieces = new_orig.rpartition('~')
            up1 = int(pieces[2])+ 1
            next_rename = pieces[0]+'~'+ str(up1)
            print(next_rename + ext)
            next_rename = rename_as_necessary(next_rename, ext)
        else:
            next_rename = new_orig + '1'
            print(next_rename + ext)
            next_rename = rename_as_necessary(next_rename, ext)
    return next_rename

##def check_duplicates(map_loc):
##    SnpFo = namedtuple('SnpFo','snp,chro,pos')
##    fo_lists = {}
##    dup_snp_list = []
##    dup_list = []
##    counter = 0
##    base, ext = os.path.splitext(map_loc)
##    with open(map_loc, mode="r") as mappy:
##        for line in mappy:
##            ls = line.split()
##            snp_fo = SnpFo(chro=ls[0],snp=ls[1],pos=ls[3])
##            chro = snp_fo.chro
##            if counter > 5000:
##                print snp_fo
##                counter = 0
##            if chro in fo_lists:
##
##                for fo in fo_lists[chro]:
##                    #if snp_fo.chro == fo.chro and
##                    if snp_fo.pos == fo.pos:
##                        print("NOTE: {0} and {1} are duplicates!"
##                              .format(snp_fo.snp, fo.snp))
##                        dup_list.append(snp_fo)
##                        dup_snp_list.append(snp_fo.snp)
##                        if fo not in dup_list:
##                            dup_list.append(fo)
##                            dup_snp_list.append(fo.snp)
##                fo_lists[chro].append(snp_fo)
##            else:
##                print('Now adding {0} as key.'.format(chro))
##                #new_list = fo_lists[chro].append(snp_fo)
##                fo_lists[chro] = list()
##                fo_lists[chro].append(snp_fo)
##                
##            counter = counter + 1
##    with open(base+'_dups.txt',mode="w") as dfile:
##        title_list = ['SNP','CHR','POS']
##        dfile.write('\t'.join(title_list) + '\n')
##        for dup in dup_list:
##            dfile.write('\t'.join(dup) +'\n')
##    return dup_snp_list
##            
##            


def usage():
    print('''
USAGE: plink_conditional.py [FLAG] OBJECT
      FLAG                  DESCRIPTION
    -m, --map               path location of binary map file 
    -h, --help              display this usage string      
    ''')
def cl_arguments(argv):
    '''
    Make use of command line arguments to assign globals.
    Args:
        argv -- sys.argv[1:], sent from main()
    '''
    global map_loc, repair_loc, dup_loc
    dup_loc = DUP_LOC
    map_loc = MAP_LOC
    repair_loc = REPAIR_LOC
    
    try: 
        opts, args = getopt.getopt(argv, "hm:d:",
                                   ["help","map=","dup=","repair="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-m","--map"):
            map_loc = arg
        elif opt in ("-d","--dup"):
            dup_loc = arg
        elif opt in ("--repair"):
            repair_loc = arg


def main(argv):
    global dup_loc, map_loc, repair_loc
    cl_arguments(sys.argv[1:])
##    with open(dup_loc, mode = "r") as dups:
##        dup_list = list()
##        for dup in dups:
##            dup_list.append(dup.strip())
    fix_map(map_loc, repair_loc)

if __name__=='__main__':
    main(sys.argv[1:])
    
