#! /ust/bin/python2.7
#! ./

'''
Created Jan 26, 2012

@author: Jessica Bonnie
'''
import os
import sys
import getopt
from operator import attrgetter

import pc_toolbox
import meta_toolbox
import fix_it

global region_loc
global table_loc
global out_file
global assoc
global build


P_BOUND = 3.5e-7
ASSOC = False

#Defines path location of allele frequency file
FREQ_LOC = '/home/jkb4y/work/data/2012Feb1/hg19/imm_controls.freq'
TABLE_LOC = '/home/jkb4y/work/data/2012Feb1/eurmeta/eurmeta_03062012_lz_hg19.txt'
REGION_LOC ='/home/jkb4y/work/data/Region_Lists/hg19/T1D_regions_hg19_04032012.txt'
OUT_FILE = '/home/jkb4y/work/results/eurmeta/hg19/InterRegionYank/inter_eurmeta.txt'

WEIGHT_MIN = 20000
BUILD = 'hg19'
blank_dict = {'chr':'--','band':'--','start':'--','end':'--','sym':'--','snp':'--',
              'p':'--','stat':'--','meta_a1':'NA','meta_a2':'NA',
              'pos':'--','note':'--','maf':'NA','control_a1':'NA',
              'control_a2':'NA','weight':''}
order_list = ['chr','band','start','end','sym','snp','pos','lz','p','z','weight',
              'meta_a1','meta_a2','maf','control_a1','control_a2','note','checkme']
def store_line(line, table_indices, assoc, annot_dict):
##    print('Now storing info.')
    if assoc:
        line_info = pc_toolbox.assoc_tuple(line, table_indices)
    else:
        line_info = meta_toolbox.meta_tuple(line, table_indices,annot_dict)
    return line_info


def read_table(table_loc, region_list, assoc, annot_dict):
    #print("Now entering read_table.")
    line1 = True
    standardizer = 1e6
    snp_list = []
    with open(table_loc, mode="r") as table:
        for line in table:
            if line1:
                if assoc:
                    table_indices = pc_toolbox.read_assoc_titles(line, .95, 'logistic')
                else:
                    table_indices = pc_toolbox.read_meta_titles(line)
                line1=False
            else:
                line_list = line.strip().split()
                snp_pos = int(line_list[table_indices['pos']])
                chro = int(line_list[table_indices['chr']])
                p_str = line_list[table_indices['p']]
                #print line_list
                p_ok = True
                if not pc_toolbox.is_number(p_str):
                    p_ok = False
                elif float(p_str)>P_BOUND:
                    p_ok=False
                in_reg = False
                for region in region_list:
                    reg_chr = int(region.chro)
                    reg_start = int(float(region.start) *standardizer)
                    reg_end = int(float(region.end) * standardizer)
##                    print("Now checking region:")
##                    print region
                    if reg_chr == chro and pc_toolbox.in_interval(snp_pos, reg_start, reg_end):
                        in_reg = True
                if p_ok and not in_reg:
                    snp_info = store_line(line, table_indices, assoc,annot_dict)
                    snp_list.append(snp_info)
    print snp_list
    sys.stdout.flush()
    tuple_list = sorted(snp_list, key=attrgetter('chro', 'p'))
    #snp_list.sort(key= lambda info : info.chro)
    print snp_list
    print tuple_list
    sys.stdout.flush()
    return tuple_list

def col_name_list(assoc,annot_dict):
    if assoc:
        title_list = ['chr','p-value','snp','bp','or','stat','ci_hi_95%','ci_lo_95%','a1']
    else:
        title_list = ['chr','p-value','im_snp','bp','lz_snp','z-score','a1','a2','weight']
    if annot_dict is not None:
        title_list.append('band')
    return title_list


def info_line(info_tuple,assoc, annot_dict):
    if assoc:
        info_list = [str(info_tuple.chro), str(info_tuple.p),
                     info_tuple.snp,str(info_tuple.pos),info_tuple.OR,
                     info_tuple.t,info_tuple.hi, info_tuple.lo,
                     info_tuple.a1]
    else:
        info_list = [str(info_tuple.chro), str(info_tuple.p),
                     info_tuple.snp,str(info_tuple.pos),info_tuple.lz,info_tuple.z,
                     info_tuple.w,info_tuple.a1, info_tuple.a2,info_tuple.band]
##    if annot_dict is not None:
##        info_list.append(annot_dict[info_tuple.band])
    text = '\t'.join(info_list) + '\n'
    return text

def write_results(info_list, out_file, assoc, annot_dict):
    if len(info_list) > 0:
        with open(out_file, mode="w") as out:
            out.write('\t'.join(col_name_list(assoc, annot_dict))+'\n')
            for info in info_list:
                text = info_line(info, assoc, annot_dict)
                out.write(text)
                print text

def usage():
    print('''
USAGE: inter_region_yank.py [FLAG] OBJECT
      FLAG                  DESCRIPTION                                 CURRENT DEFAULT
    -o, --outfile           pathname of  output file                    {0}
    -t, --table             pathname of meta file                       {1}
    -r, --region-list       pathname of region list file                {2}
        --build             genome build: hg18 or hg19                  {3}
        --assoc             (flag indicating whether table is assoc.logistic file)
    -h, --help              display this usage string



'''.format(OUT_FILE, TABLE_LOC, REGION_LOC, BUILD))


def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global table_loc, region_loc, out_file
    global assoc, fix, build
    
    table_loc = TABLE_LOC
    region_loc = REGION_LOC
    out_file = OUT_FILE
    freq_loc = FREQ_LOC
    assoc = ASSOC
    annot = None
    fix=False
    build = BUILD
    
    try: 
        opts, args = getopt.getopt(argv, "ht:r:o:f:",
                                   ["help","table=","region-list=",
                                    "out=","freq=",
                                    "assoc","build=","fix"])
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
            out_file = arg
        elif opt in ("--assoc"):
            assoc = True
        elif opt in ("--build"):
            if arg in ("hg18","hg19"):
                build = arg
            else:
                print("ERROR: build argument must be either hg18 or hg19!")
                usage()
                sys.exit(2)
        elif opt in ("--fix"):
            fix=True


def main():
    global out_file, table_loc, assoc, region_loc, build
    cl_arguments(sys.argv[1:])
    region_list = pc_toolbox.create_region_list(region_loc)
    annot_dict_loc = fix_it.locate_annot_dict(build)
    purpose = 'LOG'
    annot_dict = fix_it.build_annot_dict(purpose,annot_dict_loc)
##    annot_dict = None
##    print("Annotation Dictionary as Follows:")
##    print annot_dict
##    if annot is not None:
##        annot_dict = meta_toolbox.read_annot(annot)
##    elif fix:
##        annot_dict = meta_toolbox.read_annot_dict()
    snp_list = read_table(table_loc, region_list, assoc, annot_dict)
    print("inter_region_yank finds {0} SNPs of significance!".format(len(snp_list)))
    write_results(snp_list, out_file, assoc, annot_dict)
    
if __name__ == '__main__':
    main()
        
