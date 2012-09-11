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
from operator import itemgetter, attrgetter


def fix_map(map_loc, dup_list):
    '''
    Replaces all SNP names that do not begin with 'rs' with a name using
    chr<chromosome#>:<basepair position> format.
    Args:
        map_loc -- filepath of bim or map file
        dup_list -- list of duplicate snp names
    
    '''
    print('Now entering fix_map!')
    print('These duplicate SNPs will be left as they are:')
    print dup_list
    (basepath, ext) = os.path.splitext(map_loc)
    orig_rename = basepath + '~'
    new_loc = str(rename_as_necessary(orig_rename, ext)) + ext
    shutil.copy(map_loc, new_loc) 
    
    with open(new_loc, mode='r') as original_map:
        with open (map_loc, mode='w') as new_map:
            for orig_line in original_map:
                orig_list = orig_line.split()
                if not orig_list[1][:2]=='rs' and orig_list[1] not in dup_list:
                    new_snp = 'chr'+ orig_list[0]+':'+ orig_list[3]
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

def check_duplicates(map_loc):
    SnpFo = namedtuple('SnpFo','snp,chro,pos')
    fo_lists = {}
    dup_snp_list = list()
    dup_list = list()
    extract = list()
    counter = 0
    base, ext = os.path.splitext(map_loc)
    with open(map_loc, mode="r") as mappy:
        for line in mappy:
            ls = line.split()
            snp_fo = SnpFo(chro=ls[0],snp=ls[1],pos=ls[3])
            chro = snp_fo.chro
            if counter > 5000:
                print snp_fo
                counter = 0
            if chro in fo_lists:

                for fo in fo_lists[chro]:
                    #if snp_fo.chro == fo.chro and
                    if snp_fo.pos == fo.pos:
                        print("NOTE: {0} and {1} are duplicates!"
                              .format(snp_fo.snp, fo.snp))
                        dup_list.append(snp_fo)
                        dup_snp_list.append(snp_fo.snp)
                        if fo not in dup_list:
                            dup_list.append(fo)
                            dup_snp_list.append(fo.snp)
                fo_lists[chro].append(snp_fo)
            else:
                print('Now adding {0} as key.'.format(chro))
                #new_list = fo_lists[chro].append(snp_fo)
                fo_lists[chro] = list()
                fo_lists[chro].append(snp_fo)
                
            counter = counter + 1
    with open(base+'_dups.txt',mode="w") as dfile:
        title_list = ['SNP','CHR','POS']
        dfile.write('\t'.join(title_list) + '\n')
        for dup in dup_list:
            dfile.write('\t'.join(dup) +'\n')
    with open(base+'_extract.txt',mode="w") as efile:
        print extract
        for e in extract:
            efile.write(e +'\n')
    return dup_snp_list
            
def list_info(map_loc):
    SnpFo = namedtuple('SnpFo','snp,chro,pos')
    fo_list = list()
    extract = list()
    counter = 0
    base, ext = os.path.splitext(map_loc)
    with open(map_loc, mode="r") as mappy:
        for line in mappy:
            ls = line.split()
            snp_fo = SnpFo(chro=ls[0],snp=ls[1],pos=ls[3])
            chro = snp_fo.chro
            if counter > 5000:
                print snp_fo
                counter = 0
            fo_list.append(snp_fo)
    sort_fo = sorted(fo_list, key=attrgetter('chro', 'pos'))
    return sort_fo


def find_dups(sort_fo, extract_loc):
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
            print fo
            if not fo.snp.startswith('rs'):
                print fo
                extract.append(fo)
            if not next_fo.snp.startswith('rs'):
                print fo
                extract.append(next_fo)
        i = i + 1
        counter + 1
    print extract
    with open(extract_loc,mode="w") as efile:
        print extract
        for e in extract:
            efile.write(e.snp +'\n')
    return extract
        

def usage():
    print('''
USAGE: plink_conditional.py [FLAG] OBJECT
      FLAG                  DESCRIPTION
    -m, --map               path location of binary map file 
    -h, --help              display this usage string      
    ''')

def plink_exclude(bfile, extract_loc, script_loc):
    with open(script_loc, mode="w") as ps:
        ps.write("--bfile "+bfile+"\n")
        ps.write("--exclude "+extract_loc+"\n")
        ps.write("--out "+ bfile + "_excludeDup\n")
        ps.write("--make-bed"+"\n")
        ps.write("--noweb")

def plink_extract(bfile, extract_loc, script_loc):
    with open(script_loc, mode="w") as ps:
        ps.write("--bfile "+bfile+"\n")
        ps.write("--extract "+extract_loc+"\n")
        ps.write("--out "+ bfile + "_extractDup\n")
        ps.write("--make-bed"+"\n")
        ps.write("--noweb")

def plink_merge(exclude_b, extract_b, script_loc):
    with open(script_loc, mode="w") as ps:
        ps.write("--bfile "+exclude_b+"\n")
        ps.write("--bmerge "+extract_b+'.bed '+extract_b+'.bim '+extract_b +'.fam'+"\n")
        ps.write("--out "+ bfile + "_Dup\n")
        ps.write("--make-bed\n")
        ps.write("--noweb\n")
        ps.write("--merge-mode 6")

def plink(script_loc):
    '''runs a plink's logistic association test

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


def main(argv):
    map_loc = ''
    try: 
        opts, args = getopt.getopt(argv, "hm:",
                                   ["help","map="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-m","--map"):
            map_loc = arg
    #dup_list = check_duplicates(map_loc)
    sort_fo = list_info(map_loc)
    base, ext = os.path.splitext(map_loc)
    extract_loc = base+'_extract.txt'
    extract = find_dups(sort_fo, extract_loc)
    
    
    
    
    #fix_map(map_loc, dup_list)

if __name__=='__main__':
    main(sys.argv[1:])
    
