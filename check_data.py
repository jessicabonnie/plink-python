#! /usr/bin/python2.7
'''
Created on June 5, 2012

@author: Jessica Bonnie
'''

import sys
import os
import getopt

from collections import namedtuple
from operator import attrgetter

import fix_it
import pc_toolbox

TABLE_LOC = '/home/jkb4y/work/data/HapMap/NoDupSAMP/hg18/hapmap_removeDup~1.bim'
MANIFEST_LOC = '/home/jkb4y/work/data/Immuno_BeadChip_11419691_Manifest_B.txt'
BUILD = 'hg18'

def find_dups(sort_fo, annot_dict, repair_dict=dict(), rs_keep_list=list()):
    extract = list()
    duplicates = list()
    rs_dup = list()
    i = 0
    counter = 0
    while i < len(sort_fo)-1:
        fo = sort_fo[i]
        next_fo = sort_fo[i+1]
        if counter > 5000:
            print fo
            counter = 0
        if fo.chro == next_fo.chro and fo.pos == next_fo.pos:
            check_dup_type(annot_dict,fo,next_fo)
            if fo not in duplicates:
                duplicates.append(fo)
                duplicates.append(next_fo)
            #print("{0} and {1} are duplicates!".format(fo.snp, next_fo.snp))
            if not fo.snp.startswith('rs') and fo.snp not in repair_dict:
                #print("Adding {0} to im keep list.".format(fo.snp))
                extract.append(fo)
            elif fo.snp in rs_keep_list:
                #print("Adding {0} to im keep list.".format(fo.snp))
                extract.append(fo)
            if not next_fo.snp.startswith('rs') and next_fo.snp not in repair_dict:
                #print("Adding {0} to im keep list.".format(next_fo.snp))
                extract.append(next_fo)
            elif next_fo.snp in rs_keep_list:
                #print("Adding {0} to keep list.".format(next_fo.snp))
                extract.append(next_fo)
            if fo.snp.startswith('rs') and next_fo.snp.startswith('rs'):
##                print("FYI - {0} and {1} are rs duplicates at chrom {2}, pos {3}!"
##                      .format(fo.snp, next_fo.snp, fo.chro,fo.pos))
                pass
                if fo not in rs_dup:
                    rs_dup.append(fo)
                    rs_dup.append(next_fo)
                
        i = i + 1
        counter + 1
        
    return extract, duplicates, rs_dup

def check_dup_type(annot_dict,dup_1,dup_2):
    snp1 = dup_1.snp
    snp2 = dup_2.snp
    pair_list = [(dup_1, dup_2),(dup_2,dup_1)]
    for pair in pair_list:
        fo1 = pair[0]
        fo2 = pair[1]
        snp1 = fo1.snp
        snp2 = fo2.snp
#check if both are rs names:
        if snp1.startswith('rs') and snp2.startswith('rs'):
            if not annot_dict[snp1].rs == annot_dict[snp2].rs:
                print("SUPER WEIRD, {0} and {1} are rs duplicates at chrom {2}, pos {3}, but they don't share the same rs ids."
                      .format(snp1, snp2, fo1.chro,fo1.pos))
                print("\t{0} becomes {1}.".format(snp1,annot_dict[snp1].rs))
                print("\t{0} becomes {1}.".format(snp2,annot_dict[snp2].rs))
                return
##            print("FYI - {0} and {1} are rs duplicates at chrom {2}, pos {3}!"
##                  .format(snp1, snp2, fo1.chro,fo1.pos))
            #in the case where snp1 is not the same as it's dictionary rs name
            if not snp1 == annot_dict[snp1].rs:
                print("The rs imchip name {0} is no longer used, by rights it should now be: {1}"
                      .format(snp1, annot_dict[snp1].rs))
                if annot_dict[snp1].rs == snp2:
                    print("\tHOWEVER, this cannot be, because that is both the imchip and rs id of its duplicate.")
                elif annot_dict[snp1].rs == annot_dict[snp2].rs:
                    print("\tUNFORTUNATELY, this is also new rs id of its duplicate, {0}.".format(snp2))
            
#check if ONE is an rs name:
        elif snp1.startswith('rs') or snp2.startswith('rs'):
            if snp1.startswith('rs'):
                if annot_dict[snp1].rs == '0':
                    print("OH NOS! {0} uses an rs id that no longer exists at all!".format(snp1))
                elif not snp1 == annot_dict[snp1].rs:
                    print("The rs imchip name {0} is no longer used, by rights it should now be: {1}"
                          .format(snp1, annot_dict[snp1].rs))
                    print("\tAnd that's okay, because it's duplicate, {0}, has an imchip id.".format(snp2))
                    if not annot_dict[snp2].rs == annot_dict[snp1].rs:
                        print("\tHmm, but, it's duplicate has a different rs name listed....")
                #print("{0} is the duplicate with the rs id.".format(snp1))
                
#otherwise NEITHER is an rs name:
        else:
            print("FYI - neither {0} nor {1} has an rs-id as its original imchip id!"
                  .format(snp1, snp2))
            if annot_dict[snp1].rs == '0':
                print("\t{0} has no listed rs id.".format(snp1))
                if not annot_dict[snp2].rs == '0':
                    print("\tODDLY, {0}, on the other hand, does.".format(snp2))
                else:
                    print("\tIn fact, neither does {0}.".format(snp2))
                    return
            elif annot_dict[snp1].rs == annot_dict[snp2].rs:
                print("\tAnd they both have the same rs id, {0}.".format(annot_dict[snp1].rs))
                return
            else:
                print("WTF, {0} and {1} have different rs ids at the same position.".format(snp1,snp2))
                print("{0} has rs id {1}".format(snp1, annot_dict[snp1].rs))
                print("{0} has rs id {1}".format(snp2, annot_dict[snp2].rs))
                return
            
        


def find_dup_names(fo_list):
    duplicates = list()
    i=0
    counter = 0
    print_counter = 0
    while i < len(fo_list)-1:
        j = i+ 1
        snp = fo_list[i].snp
        fo = fo_list[i]
        print_counter = print_counter +1
        if print_counter > 1000:
            print fo
            print_counter = 0
        while j < len(fo_list):
            if fo_list[j].snp ==snp:
                print("Duplicate SNP name found: {0}".format(snp))
                if fo not in duplicates:
                    duplicates.append(fo)
                    duplicates.append(fo_list[j])
                    counter = counter + 1
            j = j + 1
        i = i + 1
    return counter, duplicates

def write_snp_list(fo_list, list_loc):
    with open(list_loc, mode="w") as ffile:
        for f in fo_list:
            ffile.write(f.snp +'\n')

def count_dups(duplicate_list, duplicate_loc):
    with open(duplicate_loc,mode="w") as dfile:
        counter = 0
        rs_counter = 0
        for d in duplicate_list:
            if d.snp.startswith('rs'):
                rs_counter = rs_counter + 1
            dfile.write(d.snp +'\n')
            counter = counter + 1
    return counter, rs_counter


def list_info(table_loc, table_type):
    SnpFo = namedtuple('SnpFo','snp,chro,pos')
    fo_list = list()
    counter = 0
    line1 = True
    with open(table_loc, mode="r") as table:
        for line in table:
            if line1:
                if table_type in ['MAP','HAPMAP']:
                    index_dict = pc_toolbox.identify_map_indices()
                elif table_type in ['META']:
                    index_dict = pc_toolbox.read_meta_titles(line)
                elif table_type == 'FAMILY':
                    index_dict = pc_toolbox.read_fam_titles(line)
                line1 = False
            else:
                ls = line.strip().split()
                snp_fo = SnpFo(chro=ls[index_dict['chr']],
                               snp=ls[index_dict['snp']],
                               pos=ls[index_dict['pos']])
                chro = snp_fo.chro
                fo_list.append(snp_fo)
    #sort_fo = sorted(fo_list, key=attrgetter('chro', 'pos'))
    return fo_list



def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global manifest_loc, table_type, build, table_loc
    manifest_loc = None
    #manifest_loc = MANIFEST_LOC
    build = BUILD
    table_loc = TABLE_LOC
    table_type = 'MAP'
    try: 
        opts, args = getopt.getopt(argv, "h",
                                   ["help","manifest=","build=","meta=",
                                    "map=","family="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("--meta"):
            table_loc = arg
            table_type = 'META'
        elif opt in ("--map"):
            table_loc = arg
            table_type = 'MAP'
        elif opt in ('--fam'):
            table_loc = arg
            table_type = 'FAMILY'
        elif opt in ("--manifest"):
            manifest_loc = arg
        elif opt in ("--build"):
            build = arg
            
def usage():
    print('''
USAGE: plink_conditional.py [FLAG] OBJECT
      FLAG                  DESCRIPTION
    --map               path location of binary map file
    --meta              path location to meta_analysis table
    -h, --help              display this usage string      
    ''')


def main(argv):
    global manifest_loc, table_loc, build, table_type, purpose
    cl_arguments(argv)
    #repair_dict = fix_it.DUPLICATE_FIX_DICT
    #rs_dup_list = fix_it.RS_KEEP
    base, ext = os.path.splitext(table_loc)
    im_dup_loc = base + '_imDups.list'
    all_dup_loc = base + '_allDups.list'
    dup_names_loc = base + '_dupNames.list'
#look up quinlan dictionary, using imchip names as key
    annot_dict_loc = fix_it.locate_annot_dict(build)
    key_type = table_type
    if not table_type == 'HAPMAP':
        key_type = 'LOG'
    annot_dict = fix_it.build_annot_dict(key_type, annot_dict_loc)

#extract SNP info from table
    print("Now gathering info from table....")
    snpFo_list = list_info(table_loc, table_type)
#check table for duplicates based on current SNP names
    print("Now searching list for duplicate snp names....")
    #rs_sort =  sorted(snpFo_list, key=attrgetter('snp'))
    dup_name_count, dup_name_list = find_dup_names(snpFo_list)
    print("FOUND {0} DUPLICATE SNP NAMES!!".format(dup_name_count))
    write_snp_list(dup_name_list, dup_names_loc)
#check table for duplicates based on current SNP positions
    print("Now sorting list for snp positions....")
    sort_fo = sorted(snpFo_list, key=attrgetter('chro', 'pos'))
    non_rs_dups, all_dups, rs_dup = find_dups(sort_fo, annot_dict)
    write_snp_list(non_rs_dups, im_dup_loc)
    
#check the duplicate list for the number of rs SNP names
    all_count, rs_count = count_dups(all_dups, all_dup_loc)
    print("\nTHERE ARE {0} TOTAL DUPLICATES --- {1} ORIGINALLY WITH RS-ID NAMES".format(all_count, rs_count))


#check for specific types of duplicates:
    #duplicates
#fix table based on current duplicate data
    #repair_dict = pc_toolbox.read_dict(repair_loc)
    repair_dict = fix_it.DUPLICATE_FIX_DICT
#recheck the duplicate values



    
if __name__=='__main__':
    main(sys.argv[1:])
    
