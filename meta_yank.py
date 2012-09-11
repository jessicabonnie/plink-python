#! /usr/bin/python2.6
#! ./
'''
Created Oct 28, 2011

@author: Jessica Bonnie
'''
import os
import re
import sys
import getopt
import pc_toolbox
import fix_it

import subprocess
import meta_toolbox
from collections import namedtuple

global region_loc
global meta_loc
global out_file
global build

TITLE_LIST = ['chr','chr_band_id','region_start(Mb)','region_end(Mb)',
              'gene_symbol','snp_name','snp_position','locuszoom_snp','imchip_name',
              'p-value','z-score', 'weight','meta_a1','meta_a2',
              'maf','control_a1','control_a2','note','check_me']


#Defines path location of allele frequency file
FREQ_LOC = None
META_LOC = None
REGION_LOC ='/home/jkb4y/work/data/Region_Lists/hg19/T1D_regions_hg19_05242012.txt.txt'
OUT_FILE = None
BFILE = None
COVAR_LOC = None
BUILD = 'hg19'

P_BOUND = 3.5e-7

WEIGHT_MIN = 10000
SEX_WEIGHT = 1400

blank_dict = {'chr':'--','ID':'--','start':'--','end':'--','sym':'--','snp':'--',
              'p':'--','stat':'--','meta_a1':'NA','meta_a2':'NA',
              'pos':'--','note':'--','maf':'NA','control_a1':'NA',
              'control_a2':'NA','weight':''}
ORDER_LIST = ['chr','ID','start','end','sym','snp','pos','lz','im','p','z','weight',
              'meta_a1','meta_a2','maf','control_a1','control_a2','note','checkme']


def usage():
    print('''
USAGE: meta_yank.py [FLAG] OBJECT
      FLAG                  DESCRIPTION                                 CURRENT DEFAULT
    -o, --outfile           pathname of  output file                    {0}
    -m, --meta              pathname of meta file                       {1}
    -f, --freq              pathname of frequency file                  {2}
    -r, --region-list       pathname of region list file                {3}
    -b, --bfile             path to plink-style binary files to use to extract OR/MAF info
        --covar             path to plink covariate file
    -h, --help              display this usage string



'''.format(OUT_FILE, META_LOC, FREQ_LOC, REGION_LOC))

def create_index_dict(table_loc, table_type):
    line1 = True
    with open(table_loc, mode = "r") as table:
        for line in table:
            if line1:
                line1 = False
                if table_type == 'family':
                    index_dict = pc_toolbox.read_fam_titles(line)
                elif table_type == 'assoc':
                    index_dict = pc_toolbox.read_assoc_titles(line, c_interval=None,plink_test='logistic')
                else:
                    index_dict = pc_toolbox.read_meta_titles(line)
    return index_dict

def read_meta(meta_loc, index_dict, gene, weight_min, table_type, annot_dict=None):
    global sex_weight, family
##    print("Weight Min is:")
##    print weight_min
    line1 = True
    standardizer = 1e6
##    if position_form == 'mb':
##        standardizer = 1e6
##    elif position_form == 'kb':
##        standardizer = 1e3
    
    old_p = 1
    old_out={}
    old_abs_z = 0
    reg_chr = gene.chro
    reg_start = int(float(gene.start) *standardizer)
    reg_end = int(float(gene.end) * standardizer)
    sym = gene.sym
    ID = gene.ID
    meta_counter = 0
    multi_counter = 1
    zero_list = []
    low_w_p = 1
    with open(meta_loc, mode="r") as meta:
        for line in meta:
            if line1:
##                if family:
##                    meta_indices = pc_toolbox.read_fam_titles(line)
##                else:
##                  meta_indices = pc_toolbox.read_meta_titles(line)
                meta_counter = 1
                line1 = False
                old_out = {'chr':'','ID':'','start':'',
                           'end':'','sym':'','snp':'','im':'',
                           'lz':'','p':'1','pos':'',
                           '|z|':'0','z':'0','note':'--',
                           'weight':'','maf':'','meta_a1':'',
                           'meta_a2':'', 'checkme':''}
                lw_blank = {'chr':'','ID':'','start':'',
                            'end':'','sym':'','snp':None,'im':None,
                            'lz':None,'p':'1','pos':'',
                            '|z|':'0','z':'0','note':'--',
                            'weight':'','maf':'','meta_a1':'',
                            'meta_a2':''}
                low_weight_out = lw_blank
                  #,'checkme':'0'}
            else:
                line_list = line.strip().split()
                snp_pos = int(line_list[index_dict['pos']])
                if int(reg_chr) == int(line_list[index_dict['chr']]
                                       ) and pc_toolbox.in_interval(snp_pos, reg_start, reg_end):
                    
                    old_p = float(old_out['p'])
                    cur_p = float(line_list[index_dict['p']])
                    lz_snp = line_list[index_dict['snp']]
                    #cur_w = float(line_list[index_dict['weight']])
                    try:
                        cur_a1 = line_list[index_dict['a1']]
                    except KeyError:
                        cur_a1 = 'NA'
                    #cur_a2 = line_list[index_dict['a2']]
                    if table_type in ['meta','family']:
                        cur_w = float(line_list[index_dict['weight']])
                        cur_a2 = line_list[index_dict['a2']]
                        cur_abs_z = abs(float(line_list[index_dict['z']]))
                        cur_z = line_list[index_dict['z']]
                        cur_snp = annot_dict[lz_snp].rs
                        cur_im = annot_dict[lz_snp].name
                        
                    else:
                        cur_w = 1
                        cur_a2 = '--'
                        cur_abs_z = float(line_list[index_dict['p']])
                        cur_z = line_list[index_dict['p']]
                        try:
                            cur_snp = annot_dict[lz_snp].rs
                            cur_im = annot_dict[lz_snp].name
                        except KeyError:
                            cur_snp = 'NA'
                            cur_im = line_list[index_dict['im']]
                    
                    #cur_snp = '--'
                    #cur_abs_z = abs(float(line_list[index_dict['z']]))
                    #cur_z = line_list[index_dict['z']]
##                    cur_snp = annot_dict[lz_snp].rs
##                    cur_im = annot_dict[lz_snp].name
                    #snp, lz_snp = meta_toolbox.correct_snp(meta_indices, line_list, annot_dict)
##                    if annot_dict is None:
##                        lz_snp=cur_snp
##                        if not cur_snp.startswith("rs"):
##                            lz_snp = 'chr'+str(reg_chr)+':'+ str(snp_pos)
##                    else:
##                        lz_snp = annot_dict[cur_snp].lz
##                    if not cur_snp.startswith('rs'):
##                        lz_snp = 'chr'+str(reg_chr)+':'+ str(snp_pos)
##                    else:
##                        lz_snp = cur_snp
                    cur_out = {'chr':gene.chro,'ID':gene.ID,'start':gene.start,
                               'end':gene.end,'sym':gene.sym,'snp':cur_snp,'im':cur_im,
                               'lz':lz_snp,'p':str(cur_p),'pos':str(snp_pos),
                               '|z|':str(cur_abs_z),'z':cur_z,'note':'--',
                               #'checkme':'0',
                               'weight':str(cur_w), 'meta_a1':cur_a1,'meta_a2':cur_a2}
                    if float(cur_out['weight']) < weight_min:
                        #print("Weight ({0}) is less than Minimum ({1}).".format(str(cur_w),weight_min))
                        cur_out['checkme']='1'
                    else:
                        cur_out['checkme']='0'
##                    if cur_w > weight_min or (reg_chr in ['23','24'] and cur_w > sex_weight):
##                        cur_out = {'chr':gene.chro,'band':gene.band,'start':gene.start,
##                                   'end':gene.end,'sym':gene.sym,'snp':cur_snp,
##                                   'lz':lz_snp,'p':str(cur_p),'pos':str(snp_pos),
##                                   '|z|':str(cur_abs_z),'z':cur_z,'note':'--','checkme':'0',
##                                   'weight':str(cur_w), 'meta_a1':cur_a1,'meta_a2':cur_a2}
##                    else:
##                        cur_out = {'chr':gene.chro,'band':gene.band,'start':gene.start,
##                                   'end':gene.end,'sym':gene.sym,'snp':cur_snp,
##                                   'lz':lz_snp,'p':str(cur_p),'pos':str(snp_pos),
##                                   '|z|':str(cur_abs_z),'z':cur_z,'note':'--','checkme':'1',
##                                   'weight':str(cur_w), 'meta_a1':cur_a1,'meta_a2':cur_a2}
                    if cur_p == old_p:
                        multi_counter = multi_counter + 1
                        ZTup = namedtuple('ZTup','snp,lz,p,z,w')
                        cur_tup = ZTup(snp=cur_out['snp'],lz=cur_out['lz'],
                                       p=cur_out['p'],z=cur_out['|z|'],w=cur_out['weight'])
                        old_tup = ZTup(snp=old_out['snp'],lz=old_out['lz'],
                                       p=old_out['p'],z=old_out['|z|'],w=old_out['weight'])
                        zero_list.append(cur_tup)
                        if old_tup not in zero_list:
                            zero_list.append(old_tup)
                        print('NOTE: Both {0} and {1} have p-values of {2}.'.format(old_out['snp'],cur_snp, cur_p))
                        print('{0} has a z-score of: {1}'.format(old_out['snp'],old_out['z']))
                        print('{0} has a z-score of: {1}'.format(cur_out['snp'],cur_out['z']))
                        if cur_out['|z|'] > old_out['|z|'] and cur_out['checkme']=='0':
                            old_p = cur_p
                            old_abs_z = cur_abs_z
                            old_out.update(cur_out)
                        print('Retaining {0} as the significant SNP based on z-score value.'.format(old_out['snp']))
                        old_out['note']='{0}SNPs(p={1})'.format(multi_counter, old_p)
                    elif cur_p < old_p:
                        if cur_out['checkme']=='0':
                            multi_counter = 1
                            #old_p = cur_p
                            old_out.update(cur_out)
                            zero_list = []
                            if not low_weight_out['p'] == '1':
                                if cur_p < float(low_weight_out['p']):
                                    low_weight_out = lw_blank
                                    low_weight_out['p']='1'
                                    low_weight_out['lz']=None
                                    print("Wiping low weight snp!")
                        elif cur_p < float(low_weight_out['p']):
                            low_weight_out.update(cur_out)
                            print('''Adding {0} to low weight SNP list due to
weight = {1}
p-value = {2}'''.format(cur_snp, low_weight_out['weight'],low_weight_out['p']))
                meta_counter = meta_counter + 1

    if old_out['chr']=='':
        print("NO SNPS FOUND IN REGION:")
        print gene
        pass
    sorted_list = sorted(zero_list, key=lambda member: member[3], reverse=True)
    #print old_out
    return old_out, sorted_list, low_weight_out['lz']

def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global meta_loc, region_loc, out_base, annot, build, fix_args
    global weight_min, sex_weight, freq_loc, covar_loc, bfile
    global family, table_type
    
    meta_loc = None
    region_loc = None
    out_base = None
    freq_loc = None
    weight_min = WEIGHT_MIN
    sex_weight = SEX_WEIGHT
    covar_loc = None
    bfile = None
    annot=None
    build = BUILD
    fix_args = list()
    family = False
    table_type = 'meta'
    
    try: 
        opts, args = getopt.getopt(argv, "hm:r:o:f:b:",
                                   ["help","meta=","region-list=",
                                    "out=","freq=","bfile=",
                                    "covar=","annot=","build=","family=","weight-min="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-f","--freq"):
            freq_loc = arg
        elif opt in ("-m","--meta"):
            meta_loc = arg
            table_type = 'meta'
            fix_args.extend(['--meta',arg])
        elif opt in ("-r","--region-list"):
            region_loc = arg
        elif opt in ("-o","--out"):
            out_base = arg
        elif opt in ("-b","--bfile"):
            bfile = arg
        elif opt in ("--family"):
            meta_loc = arg
            table_type = 'family'
            family = True
        elif opt in ("--weight-min"):
            weight_min = arg
        elif opt in ("--covar"):
            covar_loc = arg
        elif opt in ("--build"):
            build = arg
            fix_args.extend([opt,arg])
##        elif opt in ("--fix"):
##            fix = True
        elif opt in ("--annot"):
            annot = arg
       
def find_leastP_SNPs(gene_region_list, meta_loc, index_dict, out_base, freq_loc, keep_loc, weight_min,table_type,annot_dict):
    title_list = TITLE_LIST
    zero_title_list = ['snp_name','locuszoom_snp','p-value','|z-score|','weight']
    title_line = '\t'.join(title_list)
    zero_line = '\t'.join(zero_title_list)
    keep = open(keep_loc,mode="w")
    z_list = list()
    lw_list = list()
    out_file = out_base + '_yank.tbl'
    with open(out_file, mode="w") as out_text:
        out_text.write(title_line +'\n')
    with open(out_file, mode="a") as out_text:
        for gene in gene_region_list:
            print("Now searching the following region:")
            print gene
            gene_info, zero_list, lw_snp = read_meta(meta_loc, index_dict, gene, weight_min,table_type, annot_dict)#, position_form)
            keep.write(gene_info['lz']+'\n')
            if not zero_list == []:
                #zero_start, ext = os.path.splitext(out_file)
                zero_path = out_base + '_yank_'+gene_info['ID'] +'.txt'
                zeros = open(zero_path,mode='w')
                zeros.write(zero_line +'\n')
                for zero in zero_list:
                    zeros.write('\t'.join(zero)+'\n')
                zeros.close()
            freq_info = pc_toolbox.retrieve_freq(freq_loc, gene_info['snp'])
            if freq_info == None:
                freq_info = pc_toolbox.retrieve_freq(freq_loc, gene_info['lz'])
            if freq_info == None:
                gene_info.update({'maf':'NA','control_a1':'NA',
                                  'control_a2':'NA'})
            else:
                gene_info.update({'maf':freq_info[0],'control_a1':freq_info[1],
                                  'control_a2':freq_info[2]})
            for item in ORDER_LIST:
                out_text.write(str(gene_info[item]) + '\t')
            out_text.write('\n')
            if lw_snp is not None:
                lw_list.append(lw_snp)
            z_list = create_condition_list(gene_info, meta_loc, out_file, z_list)
        keep.close()
    print z_list
    print lw_list
    return z_list, lw_list
            
def create_condition_list(gene_info, meta_loc, out_file, z_list):
    head, tail = os.path.split(out_file)
    meta_head, meta_tail = os.path.split(meta_loc)
    base, ext = os.path.splitext(meta_tail)
    subfolder_path = os.path.join(head, base +'_condition_lists')
    if not os.path.exists(subfolder_path):
        os.makedirs(subfolder_path)
    filename = os.path.join(subfolder_path,gene_info['ID']+'.txt')
    clist = open(filename, mode='w')
    clist.write(gene_info['lz'])
    clist.close()
    if not gene_info['p']=='NA':
        if float(gene_info['p'])== 0:
            print("ENTERED IF STATEMENT for p==0")
            z_list.append(gene_info['lz'])
    return z_list

def plink_it(keep_loc, bfile, covar_loc, plink_out):
    put_it, tail = os.path.split(keep_loc)
    plink_cl_args = ['plink','--bfile', bfile,'--covar',covar_loc,
                     '--logistic', '--ci','.95','--noweb','--sex',
                     '--out', plink_out, '--extract', keep_loc]
    p = subprocess.Popen(plink_cl_args, bufsize = 0, executable=None,stdin=None,
                         stdout=None,stderr=None, preexec_fn=None,close_fds=False,
                         shell=False,cwd=None,env=None, universal_newlines=False,
                         startupinfo=None, creationflags=0)
    p.wait()

def repair_for_zs(z_list, lw_list, table_loc, index_dict, build):
    global family
    print ('''
LW_LIST is:''')
    print lw_list
    #fix_loc = fix_it.locate_fixed_meta(table_loc, build)
    base, ext = os.path.splitext(table_loc)
    new_loc = base + '_noZs'+ext
    line1 = True
    noZ = open(new_loc, mode = "w")
    with open(table_loc, mode="r") as tabby:
        for line in tabby:
            if line1:
##                if family:
##                    index_dict = pc_toolbox.read_fam_titles(line)
##                else:
##                    index_dict = pc_toolbox.read_meta_titles(line)
                noZ.write('\t'.join(line.strip().split())+'\n')
                line1 = False
            else:
                line_list = line.strip().split()
                snp = line_list[index_dict['snp']]
                p = line_list[index_dict['p']]
                pos = line_list[index_dict['pos']]
                if p == '0':
                    if snp in z_list:
                        line_list[index_dict['p']] = '1e-101'
                    else:
                        line_list[index_dict['p']] = '1e-100'
                elif snp in lw_list:
                    
                    print('''Changing p-value of {0} from {1} to NA due to low weight
coupled with high p-value.'''.format(snp, line_list[index_dict['p']]))
                    line_list[index_dict['p']] = 'NA'
                noZ.write('\t'.join(line_list)+'\n')
    noZ.close()
        
    
def main(argv):
    global out_base, freq_loc, meta_loc, bfile, table_type
    global covar_loc, annot, weight_min, build, fix_args
    global family
    cl_arguments(argv)
    #fix_it.main(fix_args)
    gene_region_list = pc_toolbox.create_region_list(region_loc)
    annot_dict_loc = fix_it.locate_annot_dict(build)
##    print annot_dict_loc
    annot_dict = fix_it.build_annot_dict('LOG',annot_dict_loc)
##    print("Annotation Dictionary as Follows:")
##    print annot_dict
##    gene_list = sorted(key_list,key=lambda gene:gene_region_dict[gene][0])
    #head, tail = os.path.split(out_file)
    #keep_loc = os.path.join(out_file,'_plink_keep.txt')
    keep_loc = out_base+'_plink_keep.txt'
    plink_out = out_base+'_sigs_inCC'
    #fixed_meta_loc = fix_it.locate_fixed_meta(meta_loc, build)
##    print fixed_meta_loc
    fixed_meta_loc = meta_loc
    index_dict = create_index_dict(fixed_meta_loc,table_type)
    z_list, lw_list = find_leastP_SNPs(gene_region_list, fixed_meta_loc, index_dict,
                                       out_base, freq_loc, keep_loc,
                                       weight_min, table_type, annot_dict)
    if bfile is not None and covar_loc is not None:
        plink_it(keep_loc, bfile, covar_loc, plink_out)
    repair_for_zs(z_list, lw_list, fixed_meta_loc, index_dict, build)

if __name__ == '__main__':
    main(sys.argv[1:])
    
