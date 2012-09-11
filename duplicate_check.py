#! /usr/bin/python2.7
'''
Created on Sep 14, 2011

@author: Jessica Bonnie
'''

import os
import sys
import getopt
import shutil
import subprocess
from collections import namedtuple
from operator import attrgetter
import pc_toolbox
import fix_it

REPAIR_DICT = fix_it.DUPLICATE_FIX_DICT
RS_DUP = ['rs10782171']
MESSY_LIST = ['rs3130352','1kg_6_30436336']
#FLIP_LIST = ['1kg_6_30097674', '1kg_6_30436336']
#FLIP_LIST = ['rs3130352', 'rs165255']
FLIP_LOC = '/home/jkb4y/work/data/HapMap/NoDupSAMP/DupConc'



def uniform_map(map_loc):
    '''
    Replaces all SNP names that do not begin with 'rs' with a name using
    chr<chromosome#>:<basepair position> format.
    Args:
        map_loc -- filepath of bim or map file
        dup_list -- list of duplicate snp names
    
    '''
    print('Now entering uniform_map!')
    (basepath, ext) = os.path.splitext(map_loc)
    orig_rename = basepath + '~'
    new_loc = str(rename_as_necessary(orig_rename, ext)) + ext
    shutil.copy(map_loc, new_loc) 
    
    with open(new_loc, mode='r') as original_map:
        with open (map_loc, mode='w') as new_map:
            for orig_line in original_map:
                orig_list = orig_line.split()
                
                new_snp = 'chr'+ orig_list[0]+':'+ orig_list[3]
##                if orig_list[1] not in MESSY_LIST:
##                    orig_list[1] = new_snp
                orig_list[1] = new_snp
##                if not orig_list[1][:2]=='rs' and orig_list[1] not in dup_list:
##                    new_snp = 'chr'+ orig_list[0]+':'+ orig_list[3]
##                    orig_list[1] = new_snp
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
##    dup_snp_list = list()
##    dup_list = list()
##    extract = list()
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
##    with open(base+'_extract.txt',mode="w") as efile:
##        print extract
##        for e in extract:
##            efile.write(e +'\n')
##    return dup_snp_list
            
def list_info(map_loc):
    SnpFo = namedtuple('SnpFo','snp,chro,pos')
    fo_list = list()
##    extract = list()
    counter = 0
##    base, ext = os.path.splitext(map_loc)
    with open(map_loc, mode="r") as mappy:
        for line in mappy:
            ls = line.split()
            snp_fo = SnpFo(chro=ls[0],snp=ls[1],pos=ls[3])
            chro = snp_fo.chro
##            if counter > 5000:
##                counter = 0
            fo_list.append(snp_fo)
    sort_fo = sorted(fo_list, key=attrgetter('chro', 'pos'))
    return sort_fo

def summarize_diff(merge_out):
    old_snp = 'SNP'
    counter = 1
    line1 = True
    DiffInfo = namedtuple("DiffInfo","count,snp")
    diff_list = list()
    #snp_dict = {}
    with open(merge_out+'.diff', mode="r") as diff:
        for line in diff:
            d = line.strip().split()
            new_snp = d[0]
            if not new_snp == 'SNP':
                if line1:
                    old_snp = new_snp
                    line1 = False
                elif new_snp == old_snp:
                    counter = counter + 1
                else:
                    diff_info = DiffInfo(count=counter,
                                         snp=old_snp)
                    diff_list.append(diff_info)
                    #snp_dict[old_snp]= str(counter)
                    old_snp=new_snp
                    counter = 1
        #print snp_dict
    sort_list = sorted(diff_list, reverse=True)
    print sort_list
    with open(merge_out+"_DiffSummary.txt",mode="w") as summary:
        summary.write('\t'.join(['SNP','COUNT'])+'\n')
        for diff in sort_list:
            summary.write(diff.snp + '\t'+str(diff.count)+'\n')
            

def find_dups(sort_fo, extract_loc, repair_dict):
    extract = list()
    i = 0
    counter = 0
    while i < len(sort_fo)-1:
        fo = sort_fo[i]
        next_fo = sort_fo[i+1]
        if counter > 5000:
            print fo
            counter = 0
        if fo.chro == next_fo.chro and fo.pos == next_fo.pos:
            print("{0} and {1} are duplicates!".format(fo.snp, next_fo.snp))
            if not fo.snp.startswith('rs') and fo.snp not in repair_dict:
                print("Adding {0} to extraction list.".format(fo.snp))
                extract.append(fo)
            elif fo.snp in RS_DUP:
                print("Adding {0} to extraction list.".format(fo.snp))
                extract.append(fo)
            if not next_fo.snp.startswith('rs') and next_fo.snp not in repair_dict:
                print("Adding {0} to extraction list.".format(next_fo.snp))
                extract.append(next_fo)
            elif next_fo.snp in RS_DUP:
                print("Adding {0} to extraction list.".format(next_fo.snp))
                extract.append(next_fo)
##        if fo.chro == next_fo.chro and fo.pos == next_fo.pos:
##            print("{0} and {1} are duplicates!".format(fo.snp, next_fo.snp))
##            if fo.snp.startswith('rs') and fo.snp not in RS_DUP:
##                print("Adding {0} to extraction list.".format(fo.snp))
##                extract.append(fo)
##            if next_fo.snp.startswith('rs') and next_fo.snp not in RS_DUP:
##                print("Adding {0} to extraction list.".format(next_fo.snp))
##                extract.append(next_fo)
##            if not fo.snp.startswith('rs') and fo.snp in repair_dict:
##                print("Adding {0} to extraction list.".format(fo.snp))
##                extract.append(fo)
##            if not next_fo.snp.startswith('rs')and next_fo.snp in repair_dict:
##                print("Adding {0} to extraction list.".format(next_fo.snp))
##                extract.append(next_fo)
        i = i + 1
        counter + 1
    with open(extract_loc,mode="w") as efile:
        print extract
        for e in extract:
            efile.write(e.snp +'\n')
        

def usage():
    print('''
USAGE: plink_conditional.py [FLAG] OBJECT
      FLAG                  DESCRIPTION
    -m, --map               path location of binary map file
        --flip              path location to list of SNPs to be flipped
    -h, --help              display this usage string      
    ''')

def plink_exclude(bfile, extract_loc, out, script_loc):
    with open(script_loc, mode="w") as ps:
        ps.write("--bfile "+bfile+"\n")
        ps.write("--exclude "+extract_loc+"\n")
        ps.write("--out "+ out + "\n")
        ps.write("--make-bed"+"\n")
        #ps.write("--flip "+flip_loc+"\n")
        ps.write("--noweb")

def plink_extract(bfile, extract_loc, out, script_loc,flip_loc):
    with open(script_loc, mode="w") as ps:
        ps.write("--bfile "+bfile+"\n")
        ps.write("--extract "+extract_loc+"\n")
        ps.write("--out "+ out + "\n")
        ps.write("--make-bed\n")
        ps.write("--flip "+flip_loc+"\n")
        ps.write("--noweb")

def plink_merge(exclude_b, extract_b, out, script_loc):
    with open(script_loc, mode="w") as ps:
        ps.write("--bfile "+exclude_b+"\n")
        ps.write("--bmerge "+extract_b+'.bed '+extract_b+'.bim '+extract_b +'.fam\n')
        ps.write("--out "+ out + "\n")
        #ps.write("--make-bed\n")
        ps.write("--noweb\n")
        ps.write("--merge-mode 6")

def plink(script_loc):
    '''runs a plink script 

    Keyword arguments:
    script_loc -- location of plink script

    Returns - plink subprocess
    '''
    p=subprocess.Popen(('plink','--script',script_loc),
                           bufsize = 0, executable=None,stdin=None,stdout=None,
                           stderr=None, preexec_fn=None,close_fds=False,shell=False,
                           cwd=None,env=None, universal_newlines=False,startupinfo=None,
                           creationflags=0)
    p.wait()

def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global map_loc, flip_loc, repair_loc
    map_loc = '/home/jkb4y/work/data/HapMap/NoDupSAMP/DupConc/hapmap_removeDup.bim'
    #map_loc = None
    flip_loc = FLIP_LOC
    repair_loc = None
    try: 
        opts, args = getopt.getopt(argv, "hm:",
                                   ["help","map=","flip=","repair="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-m","--map"):
            map_loc = arg
        elif opt in ("--flip"):
            flip_loc = arg
        elif opt in ("--repair"):
            repair_loc = arg

def main(argv):
    global map_loc, flip_loc, repair_loc
    cl_arguments(sys.argv[1:])
##    map_loc = ''
##    try: 
##        opts, args = getopt.getopt(argv, "hm:",
##                                   ["help","map="])
##    except getopt.GetoptError:
##        usage()
##        sys.exit(2)
##    for opt, arg in opts:
##        if opt in ("-h","--help"):
##            usage()
##            sys.exit()
##        elif opt in ("-m","--map"):
##            map_loc = arg
##        elif opt in ("--flip"):
##            flip_loc = arg

    repair_dict = pc_toolbox.read_dict(repair_loc)
    
    sort_fo = list_info(map_loc)
    base, ext = os.path.splitext(map_loc)
    extract_loc = base+'_extractList.txt'
    exclude_bfile = base + "_excludeDup"
    extract_bfile = base + "_extractDup"
    find_dups(sort_fo, extract_loc, repair_dict)
    
    #write and use the script to extract the non-rs duplicates
    extract_script = base + '_extractScript.txt'
    plink_extract(base, extract_loc, extract_bfile,extract_script,flip_loc)
    plink(extract_script)
    
    #write and use the script to exclude the non-rs duplicates
    exclude_script = base + '_excludeScript.txt'
    plink_exclude(base, extract_loc, exclude_bfile,exclude_script)
    plink(exclude_script)
    
    #convert the map files from both results to use only chr:position form
    uniform_map(extract_bfile+'.bim')
    uniform_map(exclude_bfile+'.bim')
    #run concordance check
    merge_script = base + '_mergeScript.txt'
    merge_out = base + '_DupMerge'
    plink_merge(exclude_bfile, extract_bfile,merge_out, merge_script)
    plink(merge_script)

    #summarize the results of the diff file
    summarize_diff(merge_out)
    #if the concordance check comes out perfectly, then keep the _exclude~.bim as final
    
    
    #fix_map(map_loc, dup_list)

if __name__=='__main__':
    main(sys.argv[1:])
    
