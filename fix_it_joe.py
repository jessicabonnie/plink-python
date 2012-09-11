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
import fix_it

import meta_toolbox
import pc_toolbox
import look_up_rs_joe

BUILD = 'hg18'



def fix_joe_map(table_loc, annot_dict,build):
    counter = 0
    base, ext = os.path.splitext(table_loc)
    new_loc = base + '_06082012'+ext
    with open(table_loc, mode="r") as table:
        output = open(new_loc, mode="w")
        index_dict = {'chr':0,'pos':3,'snp':1}
        error_list = list()
        for line in table:
            line_list = line.strip().split()
            snp = line_list[index_dict['snp']]
            if counter > 100:
                print line_list
                counter = 0
            if build == 'hg18':
                try:
                    if not line_list[index_dict['chr']] == annot_dict[snp].hg18_chr:
                        print("{0} has mismatched chr info!".format(snp))
                    if not line_list[index_dict['pos']] == annot_dict[snp].hg18_pos:
                        print("{0} has mismatched pos info!".format(snp))
                    if not annot_dict[snp].rs == '0':
                        line_list[index_dict['snp']]= annot_dict[snp].rs
                except KeyError:
                    error_list.append(list(line_list))
                    print("{0} is missing from the dictionary!".format(line_list[index_dict['snp']]))
            if build == 'hg19':
                try:
                    if not line_list[index_dict['chr']] == annot_dict[snp].hg18_chr:
                        print("{0} has mismatched chr info!".format(snp))
                    if not line_list[index_dict['pos']] == annot_dict[snp].hg18_pos:
                        print("{0} has mismatched pos info!".format(snp))
                    if not annot_dict[snp].rs == '0':
                        line_list[index_dict['snp']]= annot_dict[snp].rs
                except KeyError:
                    error_list.append(list(line_list))
                    print("{0} is missing from the dictionary!".format(line_list[index_dict['snp']]))
            output.write('\t'.join(line_list)+'\n')
            counter = counter + 1
    output.close()
    base, ext = os.path.splitext(table_loc)
    error_loc = base + '_NAMEERRORS.txt'
    efile = open(error_loc,mode="w")
    efile.write('\t'.join(['CHR','RS','cM','POS','A1','A2'])+'\n')
    for error in error_list:
        efile.write('\t'.join(error)+'\n')
    efile.close()
    return True
                





def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global annot_loc, table_loc, meta, build, hapmap, family
    table_loc = look_up_rs.MAP_FILE
    build = BUILD
    hapmap = False
    family = False
    try: 
        opts, args = getopt.getopt(argv, "h",
                                   ["help","annot=","map=","meta=","build=",
                                    "hapmap","family="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("--map"):
            table_loc = arg
        elif opt in ("--meta"):
            table_loc = arg
            meta = True
        elif opt in ("--family"):
            table_loc = arg
            family = True
        elif opt in ("--hapmap"):
            hapmap = True
        elif opt in ("--annot"):
            annot_loc = None
        elif opt in ("--build"):
            build = arg
        
def usage():
    print('''
USAGE: fix_it.py [FLAG] OBJECT
      FLAG                  DESCRIPTION                         CURRENT DEFAULT
        --meta              filepath of meta table                      
        --map               filepath of plink map file                  
        --annot             filepath of annotation file         {0}                
        --build             genome build ( hg18 or hg19 )       {1}
    -h, --help              display this usage string



'''.format(BUILD, ANNOT_LOC))



def main(argv):
    global table_loc, build
    cl_arguments(argv)
    build = 'hg18'

    dict_loc = fix_it.locate_annot_dict(build)
    annot_dict = fix_it.build_annot_dict('MAP', dict_loc)
    tf = fix_joe_map(table_loc,annot_dict,build)
    look_up_rs_joe.main()

    
if __name__=='__main__':
    main(sys.argv[1:])
