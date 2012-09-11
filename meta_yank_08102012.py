#! /usr/bin/python2.6
#! ./
'''
Created Aug 10, 2012

@author: Jessica Bonnie
'''

#IMPORTS
import os
import sys
import getopt
import subprocess
import itertools
from collections import namedtuple

import pc_toolbox
import fix_it


#DEFAULTS
REGION_LOC ='/home/jkb4y/work/data/Region_Lists/hg19/T1D_regions_hg19_05242012.txt'
BUILD = 'hg19'
P_BOUND = 3.5e-7
WEIGHT_MIN = 10000
#
META_TITLE_LIST = ['chr','chr_band_id','region_start(Mb)','region_end(Mb)',
                   'chr_band','snp_name','snp_position','locuszoom_snp',
                   'imchip_name','p-value','z-score', 'weight','meta_a1',
                   'meta_a2','maf','control_a1','control_a2','note']
META_ORDER_LIST = ['chr','ID','start','end','band','snp','pos','lz','im',
                   'p','z','weight','meta_a1','meta_a2','maf',
                   'control_a1','control_a2','note']
PERM_ORDER_LIST = ['ID','chr','snp_ap','pos','p','gene','note','np','maf',
                   'control_a1','control_a2','snp_jb','im']
PERM_TITLE_LIST = ['chr_band_id','chr','SNP_AP','snp_position','p-value','GENE','NOTE',
                   'NP','MAF_UK','CONTROL_A1','CONTROL_A2','locuszoom_snp','imchip_name']

ASSOC_TITLE_LIST = ['chr','chr_band_id','region_start(Mb)','region_end(Mb)',
                    'chr_band','snp_name','snp_position','locuszoom_snp',
                    'imchip_name','p-value','stat','OR','a1','ci_lo','ci_hi',
                    'maf','control_a1','control_a2','note']
ASSOC_ORDER_LIST = ['chr','ID','start','end','band','snp','pos','lz','im','p',
                    't','or','a1','hi','lo','maf','control_a1','control_a2','note']

AA_ASSOC_TITLE_LIST = ['chr','chr_band_id','region_start(Mb)','region_end(Mb)',
                       'chr_band','snp_name','snp_position','locuszoom_snp',
                       'imchip_name','p-value','STAT','OR','A1','ci_hi','ci_lo',
                       'AA_maf','AA_control_a1','AA_control_a2',
                       'UK_maf','UK_control_a1','UK_control_a2','note']
AA_ASSOC_ORDER_LIST = ['chr','ID','start','end','band','snp','pos','lz','im','p',
                       't','or','a1','hi','lo',
                       'AA_maf','AA_control_a1','AA_control_a2',
                       'UK_maf','UK_control_a1','UK_control_a2','note']


EQTL_ORDER_LIST = ['ID','chr','snp_ap','pos','p','t','gene','note','maf',
                   'control_a1','control_a2','snp_jb','im']

EQTL_TITLE_LIST = ['CHR_BAND_ID','CHR','SNP_AP','BP','P','STAT','GENE','NOTE',
                   'NP','MAF','CONTROL_A1','CONTROL_A2','SNP_JB','SNP_IM']

ZTup = namedtuple('ZTup','snp,p,abs_z,w')




#read the line and compare it's values with those in the argument
def compare_line(line, old_line, index_dict, table_type, multicounter, lw_list,weight_min = 0):
    line_list = line.strip().split()
    current_sig = old_line.strip().split()
    checkme = '0'
    cur_snp = line_list[index_dict['snp']]
    old_p = float(current_sig[index_dict['p']])
    cur_p = float(line_list[index_dict['p']])
    if table_type in ['meta','family']:
    #only meta and family tables have z-scores or weights to compare
        cur_z = line_list[index_dict['z']]
        cur_abs_z = abs(float(cur_z))
        old_comp = abs(float(current_sig[index_dict['z']]))
        second_comp = abs(float(cur_z))
        cur_w = float(line_list[index_dict['weight']])
        if cur_w < weight_min:
            checkme ='1'
            lw_list.append(cur_snp)
    elif table_type in ['assoc','eqtl']:
        try:
            old_comp = abs(float(current_sig[index_dict['t']]))
            cur_t = line_list[index_dict['t']]
            second_comp = abs(float(cur_t))
        except KeyError:
            old_comp = abs(float(current_sig[index_dict['p']]))
            cur_t = line_list[index_dict['p']]
            second_comp = abs(float(cur_t))
    else:
        old_comp = 1
        second_comp = 1
    if cur_p == old_p:
        print('NOTE: Both {0} and {1} have p-values of {2}.'.format(current_sig[index_dict['snp']],
                                                                    cur_snp,
                                                                    cur_p))
        print('{0} has a z-score/stat of: {1}'.format(current_sig[index_dict['snp']],old_comp))
        print('{0} has a z-score/stat of: {1}'.format(cur_snp,second_comp))
        multicounter = multicounter + 1
        if second_comp > old_comp:
            print('Retaining {0} as the significant SNP based on z-score/stat value.'
                  .format(cur_snp))
            return True, multicounter, lw_list
        else:
            print('Retaining {0} as the significant SNP based on z-score/stat value.'
                  .format(current_sig[index_dict['snp']]))
            return False, multicounter, lw_list
    elif cur_p < old_p:
        multicounter = 1
        return True, multicounter, lw_list
    else:
        return False, multicounter, lw_list

def assoc_line_dict(line_list, assoc_indices):
    assoc_dict = dict()
    assoc_dict['t']=line_list[assoc_indices['t']]
    return assoc_dict
    

def make_info_dict(line_list, region, index_dict, table_type, log_annot_dict, pop, freq_loc, aa_freq_loc=None):
    info_dict = dict()
    snp_key = line_list[index_dict['snp']]
    if snp_key == 'NA':
        print info_dict
    info_dict['snp'] = snp_key
    info_dict['p'] = line_list[index_dict['p']]
    info_dict['pos'] = line_list[index_dict['pos']]
    print(info_dict)
    if table_type in ['family','meta']:
        type_dict = meta_line_dict(line_list, index_dict)
    elif table_type in ['eqtl']:
        type_dict = eqtl_line_dict(line_list, index_dict)
    elif table_type in ['assoc']:
        type_dict = assoc_line_dict(line_list, index_dict)
    elif table_type in ['perm']:
        type_dict = perm_line_dict(line_list, index_dict)
    print("5")
    info_dict.update(type_dict)
    print("6")
    print info_dict
    info_dict = update_with_region(region, info_dict)
    #info_dict = update_with_freq(freq_loc, info_dict)
    if table_type in ['eqtl','perm']:
        info_dict = update_eqtl_snp_names(line_list, index_dict, log_annot_dict, info_dict)
    else:
        info_dict = update_with_snp_names(log_annot_dict, snp_key, info_dict)
    info_dict = update_with_freq(freq_loc, info_dict, pop, aa_freq_loc)
    #info_dict = update_with_freq(freq_loc, info_dict)
    return info_dict

def eqtl_line_dict(line_list, eqtl_indices):
    eqtl_dict = dict()
    eqtl_dict['t']=line_list[eqtl_indices['t']]
    eqtl_dict['gene']=line_list[eqtl_indices['gene']]
    eqtl_dict['a1']=line_list[eqtl_indices['a1']]
    eqtl_dict['or']=line_list[eqtl_indices['or']]
    return eqtl_dict

def perm_line_dict(line_list, perm_indices):
    perm_dict = dict()
    perm_dict['gene']=line_list[perm_indices['gene']]
    perm_dict['np']=line_list[perm_indices['np']]
    return perm_dict

def dict_helper(line_list, index):
    '''Deals with situations in which there is no index, because there
        is no column containing the specific information.
    '''
    item = ''
    if index is not None:
        item = line_list[int(index)]
    sys.stdout.flush()
    return item

def assoc_line_dict(line_list, assoc_indices):
    assoc_dict = dict()
    assoc_dict['lo']=dict_helper(line_list,assoc_indices['lo'])
    assoc_dict['hi']=dict_helper(line_list,assoc_indices['hi'])
    assoc_dict['t']=dict_helper(line_list,assoc_indices['t'])
    assoc_dict['or']=dict_helper(line_list,assoc_indices['or'])
    assoc_dict['a1']=line_list[assoc_indices['a1']]
    return assoc_dict
    

def meta_line_dict(line_list,  meta_indices):
    meta_dict = dict()
##    lz_snp = line_list[meta_indices['snp']]
##    meta_dict['p'] = line_list[meta_indices['p']]
    meta_dict['z']=line_list[meta_indices['z']]
    meta_dict['weight'] = line_list[meta_indices['weight']]
##    meta_dict['pos'] = line_list[meta_indices['pos']]
    meta_dict['meta_a1'] = line_list[meta_indices['a1']]
    meta_dict['meta_a2'] = line_list[meta_indices['a2']]
    
    return meta_dict

def update_with_region(region, snp_info):
    snp_info.update({'chr':region.chro,'start':region.start,
                     'end':region.end, 'ID':region.ID, 'band':region.title})
    return snp_info

def update_eqtl_snp_names(line_list, index_dict, annot_dict, snp_info):
    print("1")
    print line_list[index_dict['snp']]
    if line_list[index_dict['snp']] == 'NA':
        print("2")
        snp_info.update({'snp':'NA',
                         'im':'NA',
                         'lz':'NA',
                         'snp_jb':'NA',
                         'snp_ap':'NA'})
    else:
        im_name = line_list[index_dict['im']]
        jb_name = line_list[index_dict['snp']]
        lz_name = line_list[index_dict['lz']]
        ap_name = line_list[index_dict['snp_ap']]
        try:
            rs_name = annot_dict[jb_name].rs
        except KeyError:
            rs_name = 'NA'
        snp_info.update({'snp':rs_name,
                         'snp_ap':ap_name,
                         'im':im_name,
                         'snp_jb':jb_name,
                         'lz':jb_name})
    print snp_info
    return snp_info

def update_with_snp_names(annot_dict, snp_key, snp_info):
    print snp_info
    rs_name = annot_dict[snp_key].rs
    im_name = annot_dict[snp_key].name
    lz_name = annot_dict[snp_key].lz
    
    snp_info.update({'snp':rs_name,
                     'im':im_name,
                     'lz':lz_name})
    return snp_info

def update_with_freq(freq_loc, snp_info, pop, aa_freq_loc=None):
    print snp_info
    print ("3")
    if pop == 'UK' and aa_freq_loc == None:
        if snp_info['snp'] == 'NA':
            snp_info.update({'maf':'NA','control_a1':'NA','control_a2':'NA'})
            return snp_info
        freq_info = pc_toolbox.retrieve_freq(freq_loc, snp_info['snp'])
        if freq_info == None:
            freq_info = pc_toolbox.retrieve_freq(freq_loc, snp_info['lz'])
        if freq_info == None:
            snp_info.update({'maf':'NA','control_a1':'NA','control_a2':'NA'})
        else:
            snp_info.update({'maf':freq_info[0],'control_a1':freq_info[1],
                             'control_a2':freq_info[2]})
        return snp_info
    else:
        snp_info = update_with_aa_freq(aa_freq_loc, freq_loc, snp_info)
        return snp_info

def update_with_aa_freq(aa_freq_loc,uk_freq_loc, snp_info):
    print snp_info
    print ("3")
    if snp_info['snp'] == 'NA':
        snp_info.update({'AA_maf':'NA','AA_control_a1':'NA','AA_control_a2':'NA'})
        return snp_info
    aa_freq_info = pc_toolbox.retrieve_freq(aa_freq_loc, snp_info['snp'])
    if aa_freq_info == None:
        aa_freq_info = pc_toolbox.retrieve_freq(aa_freq_loc, snp_info['lz'])
    if aa_freq_info == None:
        snp_info.update({'AA_maf':'NA','AA_control_a1':'NA','AA_control_a2':'NA'})
    else:
        snp_info.update({'AA_maf':aa_freq_info[0],'AA_control_a1':aa_freq_info[1],
                          'AA_control_a2':aa_freq_info[2]})
    uk_freq_info = pc_toolbox.retrieve_freq(uk_freq_loc, snp_info['snp'])
    if uk_freq_info == None:
        uk_freq_info = pc_toolbox.retrieve_freq(uk_freq_loc, snp_info['lz'])
    if uk_freq_info == None:
        snp_info.update({'UK_maf':'NA','UK_control_a1':'NA','UK_control_a2':'NA'})
    else:
        snp_info.update({'UK_maf':uk_freq_info[0],'UK_control_a1':uk_freq_info[1],
                          'UK_control_a2':uk_freq_info[2]})
    print snp_info
    return snp_info


def standardize_region(region):
    standardizer = 1e6
    reg_chr = region.chro
    reg_start = int(float(region.start) *standardizer)
    reg_end = int(float(region.end) * standardizer)
    return (reg_chr, reg_start, reg_end)
    
    
#read the line_list and determine if it is in the interval
def snp_in_interval(region_tuple, line_list, index_dict):
    reg_chr,reg_start,reg_end = region_tuple
    snp_pos = int(line_list[index_dict['pos']])
    if not int(reg_chr) == int(line_list[index_dict['chr']]):
        return False
    elif pc_toolbox.in_interval(snp_pos, reg_start, reg_end):
        return True
    else:
        return False
def make_base_line(index_dict, table_type):
    length = len(index_dict.keys())
    base_list = list(itertools.repeat("--", length+1))
    base_list[index_dict['snp']] = 'NA'
    base_list[index_dict['pos']] = '0'
    base_list[index_dict['p']] = '1'
    if table_type in ['meta','family']:
        base_list[index_dict['z']] = '0'
    elif table_type in ['assoc','eqtl']:
        base_list[index_dict['t']] = '0'
    if table_type in ['eqtl','perm']:
        base_list[index_dict['gene']] = 'NA'
    base_line = '\t'.join(base_list)
    return base_line


def meta_zlist(line, index_dict, zlist):
    line_list = line.strip().split()
    snp = line_list[index_dict['snp']]
    p = float(line_list[index_dict['p']])
    abs_z = abs(float(line_list[index_dict['z']]))
    w = float(line_list[index_dict['weight']])
    if p == 0:
        ztup = ZTup(snp=snp,p=p,abs_z=abs_z,w=w)
        zlist.append(ztup)
    return zlist


def sort_zlist(zlist):
    sorted_list = sorted(zlist, key=lambda member: member[2], reverse=True)
    print sorted_list
    return sorted_list
    
            

def read_table(table_loc, index_dict, region, weight_min, table_type, annot_dict, zlist):
    line1 = True
    old_p = 1
    old_abs_z = 0
    lw_list = list()
    region_tuple = standardize_region(region)
    print region_tuple
    (reg_chr, reg_start, reg_end) = region_tuple
    band = region.band
    ID = region.ID
    multicounter = 1
    empty_line = make_base_line(index_dict, table_type)
    with open(table_loc, mode="r") as table:
        for line in table:
            if line1:
                base_line = empty_line
                line1 = False
            else:
                line_list = line.strip().split()
                if snp_in_interval(region_tuple, line_list,index_dict):
                    significant, multicounter, lw_list = compare_line(line, base_line, index_dict,
                                                             table_type, multicounter,
                                                             weight_min)
                    if significant:
                        base_line = line
    if table_type in ['meta','family']:
        zlist = meta_zlist(base_line,index_dict, zlist)
##                        snp_info.clear()
##                        snp_info = make_info_dict(line_list, region, index_dict, table_type, annot_dict)
    return base_line, zlist, multicounter


def create_annot_dict(build, purpose):
    annot_dict_loc = fix_it.locate_annot_dict(build)
    annot_dict = fix_it.build_annot_dict(purpose,annot_dict_loc)
    return annot_dict

##def create_index_dict(table_loc, table_type):
##    line1 = True
##    with open(table_loc, mode = "r") as table:
##        for line in table:
##            if line1:
##                line1 = False
##                if table_type == 'meta':
##                    index_dict = pc_toolbox.read_meta_titles(line)
##                elif table_type == 'family':
##                    index_dict = pc_toolbox.read_fam_titles(line)
##                elif table_type == 'assoc':
##                    index_dict = pc_toolbox.read_assoc_titles(line, c_interval=.95,plink_test='logistic')
##                elif table_type == 'eqtl':
##                    index_dict = pc_toolbox.read_eqtl_titles(line)
##                elif table_type == 'perm':
##                    index_dict = pc_toolbox.read_perm_titles(line)
##                return index_dict, line

def document_result(sig_line, region, index_dict, table_type,annot_dict, multicounter):
    sig_line_list = sig_line.strip().split()
    snp_info = make_info_dict(sig_line_list, region, index_dict, table_type,
                              annot_dict,pop, freq_loc, aa_freq_loc)
    print("4")
    print snp_info
    snp_info['note'] = '--'
    if multicounter > 1:
        snp_info['note'] = '{0}SNPs(p={1})'.format(multicounter, snp_info['p'])
    if snp_info['p']=='1':
        print("NO SNPS FOUND IN REGION:")
        print region
    print snp_info
    return snp_info


def repair_meta_for_lz(zlist, table_loc, index_dict, weight_min):
    print ('''
Z_LIST is:''')
    print zlist
    if len(zlist) == 0:
        return
    zeros = list()
    for z in zlist:
        zeros.append(z.snp)
    base, ext = os.path.splitext(table_loc)
    new_loc = base + '_noZs'+ext
    line1 = True
    noZ = open(new_loc, mode = "w")
    with open(table_loc, mode="r") as tabby:
        for line in tabby:
            line_list = line.strip().split()
            if line1:
                noZ.write('\t'.join(line_list)+'\n')
                line1 = False
            else:
                snp = line_list[index_dict['snp']]
                p = line_list[index_dict['p']]
                if float(p) == 0:
                    if snp in zeros:
                        line_list[index_dict['p']] = '1e-101'
                    else:
                        line_list[index_dict['p']] = '1e-100'
                if float(line_list[index_dict['weight']]) < weight_min:
                    line_list[index_dict['p']] = 'NA'
                noZ.write('\t'.join(line_list)+'\n')
    noZ.close()


def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global table_loc, region_loc, out_base, annot, build
    global weight_min, freq_loc, covar_loc, bfile
    global table_type, pop, aa_freq_loc
    
    table_loc = None
    region_loc = REGION_LOC
    out_base = None
    freq_loc = None
    weight_min = WEIGHT_MIN
    covar_loc = None
    bfile = None
    annot=None
    build = BUILD
    table_type = 'meta'
    pop = "UK"
    aa_freq_loc = None
    
    try: 
        opts, args = getopt.getopt(argv, "hm:r:o:f:b:",
                                   ["help","meta=","region-list=",
                                    "out=","freq=","bfile=",
                                    "covar=","annot=","build=","family=",
                                    "weight-min=","eqtl=","perm=","assoc=",
                                    "pop=","freq-aa="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-f","--freq"):
            freq_loc = arg
        elif opt in ("--freq-aa"):
            aa_freq_loc = arg
        elif opt in ("-m","--meta"):
            table_loc = arg
            table_type = 'meta'
        elif opt in ("-r","--region-list"):
            region_loc = arg
        elif opt in ("-o","--out"):
            out_base = arg
        elif opt in ("-b","--bfile"):
            bfile = arg
        elif opt in ("--family"):
            table_loc = arg
            table_type = 'family'
        elif opt in ("--eqtl"):
            table_loc = arg
            table_type = 'eqtl'
        elif opt in ("--assoc"):
            table_loc = arg
            table_type = 'assoc'
        elif opt in ("--perm"):
            table_loc = arg
            table_type = 'perm'
        elif opt in ("--weight-min"):
            weight_min = arg
        elif opt in ("--covar"):
            covar_loc = arg
        elif opt in ("--build"):
            build = arg
        elif opt in ("--pop"):
            pop = arg
        elif opt in ("--annot"):
            annot = arg
            
def create_output_title_line(table_type, input_title_line):
    line_list = input_title_line.strip().split()
    if table_type in ['meta','family']:
        out_line = '\t'.join(META_TITLE_LIST)
    if table_type == 'perm':
        out_line = '\t'.join(PERM_TITLE_LIST)
    if table_type == 'eqtl':
        line_list.extend(['note','chr_band','band_id'])
        out_line = '\t'.join(line_list)
    if table_type == 'assoc':
        if pop == 'AA':
            out_line = '\t'.join(AA_ASSOC_TITLE_LIST)
        if pop == 'UK':
            out_line = '\t'.join(ASSOC_TITLE_LIST)
    return out_line
        

def main(argv):
    global out_base, freq_loc, table_loc, bfile, table_type
    global covar_loc, annot, weight_min, build, region_loc, pop, aa_freq_loc
    cl_arguments(argv)
    zlist = list()
    index_dict, input_title_line = pc_toolbox.create_index_dict(table_loc, table_type)
    region_list = pc_toolbox.create_region_list(region_loc)
    annot_dict = create_annot_dict(build, 'LOG')
    out_file = out_base + '_yank.tbl'
    output_title_line = create_output_title_line(table_type, input_title_line)
    with open(out_file, mode="w") as out_text:
##        title_line = '\t'.join(META_TITLE_LIST)
        out_text.write(output_title_line +'\n')
        for region in region_list:
            print zlist
            sig_line, zlist, multicounter = read_table(table_loc, index_dict, region, weight_min, table_type, annot_dict, zlist)
            info = document_result(sig_line, region, index_dict, table_type,annot_dict, multicounter)
            if table_type == 'assoc':
                if pop == 'AA':
                    for item in AA_ASSOC_ORDER_LIST:
                        out_text.write(str(info[item]) + '\t')
                if pop == 'UK':
                    for item in ASSOC_ORDER_LIST:
                        out_text.write(str(info[item]) + '\t')
                
##            if table_type in ['perm']:
##                out_list = sig_line.strip().split()
##                out_list.append(info['note'])
##                out_list.append(info['band'])
##                out_list.append(info['ID'])
##                out_text.write('\t'.join(out_list))
            if table_type in ['meta','family']:
                for item in META_ORDER_LIST:
                    out_text.write(str(info[item]) + '\t')
            if table_type =='perm':
                for item in PERM_ORDER_LIST:
                    out_text.write(str(info[item]) + '\t')
            out_text.write('\n')
    if table_type in ['meta','family']:
        repair_meta_for_lz(zlist, table_loc, index_dict, weight_min)
    
if __name__ == '__main__':
    main(sys.argv[1:])
