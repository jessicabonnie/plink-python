#! /usr/bin/python2.7
#! ./
'''
Created on Dec 14, 2011

@author Jessica Bonnie
'''
from collections import namedtuple
import sys
import os

import fix_it

RegInfo = namedtuple('RegInfo', 'chro,start,end,sym,band,title,ID')

REG_CHR_TITLE = 'gene_chr'
START_TITLE = 'region_start'
END_TITLE = 'region_end'
GENE_TITLE = 'gene_symbol'
BAND_TITLE = 'Chr_band'
TITLE_TITLE = 'band_title'
ID_TITLE='region_id'


META_CHR_TITLE = 'CHR'
SNP_TITLE = 'MarkerName'
P_TITLE = 'P-value'
POSITION_TITLE = 'POS'
Z_TITLE = 'Zscore'
WEIGHT_TITLE = 'Weight'
A1_TITLE = 'Allele1'
A2_TITLE = 'Allele2'
DIRECTION_TITLE = 'Direction'
FLAG_TITLE = 'eurlowP_SNP'

FAM_CHR_TITLE = 'CHR'
FAM_SNP_TITLE = 'SNP'
FAM_P_TITLE = 'PVALUE'
FAM_POSITION_TITLE = 'POS'
FAM_Z_TITLE = 'Z_GDT'
FAM_WEIGHT_TITLE = 'Pair'
FAM_A1_TITLE = 'AL1'
FAM_A2_TITLE = 'AL2'

YANK_LZ_TITLE = 'locuszoom_snp'
YANK_POS_TITLE = 'snp_position'
YANK_IM_TITLE = 'imchip_name'
YANK_REG_START_TITLE = 'region_start(Mb)'
YANK_REG_END_TITLE = 'region_end(Mb)'
YANK_BAND_ID_TITLE = 'chr_band_id'
YANK_CHR_TITLE = 'chr'
YANK_P_TITLE = 'p-value'

JB_SNP_TITLE = 'SNP_LZ_JB'
IM_SNP_TITLE = 'SNP_IM'
AP_SNP_TITLE = 'SNP_AP'


def read_eqtl_titles(title_line):
    title_list = title_line.strip().split()
    a1_index = None
    ci_hi_index = None
    ci_lo_index = None
    or_index = None
    stat_index = None
##    lz_choice = JB_SNP_TITLE
##    print lz_choice
##    print JB_SNP_TITLE
    p_index = title_list.index('P')
    a1_index = title_list.index('A1')
    or_index = title_list.index('BETA')
    stat_index = title_list.index('STAT')
    gene_index = title_list.index('GENE')
    index_dict = {'p':p_index,
                  'chr':title_list.index('CHR'),
                  'pos':title_list.index('BP'),
                  'snp':title_list.index(JB_SNP_TITLE),
                  'snp_ap':title_list.index(AP_SNP_TITLE),
                  'im':title_list.index(IM_SNP_TITLE),
                  'or':or_index,
                  't':stat_index,
                  'hi':ci_hi_index,
                  'lo':ci_lo_index,
                  'a1':a1_index,
                  'pvalcol':'P',
                  'markercol':lz_choice,
                  'gene':gene_index}
    return index_dict




def create_index_dict(table_loc, table_type, c_interval=.95):
    line1 = True
    with open(table_loc, mode = "r") as table:
        for line in table:
            if line1:
                line1 = False
                if table_type == 'meta':
                    index_dict = read_meta_titles(line)
                elif table_type == 'family':
                    index_dict = read_fam_titles(line)
                elif table_type == 'assoc':
                    index_dict = read_assoc_titles(line, c_interval,plink_test='logistic')
                elif table_type == 'eqtl':
                    index_dict = read_eqtl_titles(line)
                elif table_type == 'perm':
                    index_dict = read_perm_titles(line)
                return index_dict, line

def read_perm_titles(title_line):
    title_list = title_line.strip().split()
    a1_index = None
##    ci_hi_index = None
##    ci_lo_index = None
##    or_index = None
##    stat_index = None
    lz_choice = JB_SNP_TITLE
    p_index =  title_list.index('EMP1')
    gene_index = title_list.index('GENE')
    np_index = title_list.index('NP')
    index_dict = {'p':p_index,
                  'chr':title_list.index('CHR'),
                  'pos':title_list.index('BP'),
                  'snp':title_list.index(JB_SNP_TITLE),
                  'lz':title_list.index(lz_choice),
                  'snp_ap':title_list.index(AP_SNP_TITLE),
                  'im':title_list.index(IM_SNP_TITLE),
##                  'or':or_index,
##                  't':stat_index,
##                  'hi':ci_hi_index,
##                  'lo':ci_lo_index,
                  'a1':a1_index,
                  'np':np_index,
                  'pvalcol':'EMP1',
                  'markercol':lz_choice,
                  'gene':gene_index}
    return index_dict


def read_fam_titles(line):
    print line
    title_list = line.strip().split()
    meta_chr_col = title_list.index(FAM_CHR_TITLE)
    p_col = title_list.index(FAM_P_TITLE)
    snp_col = title_list.index(FAM_SNP_TITLE)
    position_col = title_list.index(FAM_POSITION_TITLE)
    z_col = title_list.index(FAM_Z_TITLE)
    weight = title_list.index(FAM_WEIGHT_TITLE)
    meta_a = title_list.index(FAM_A1_TITLE)
    a2 = title_list.index(FAM_A2_TITLE)
    fam_indices = {'chr':meta_chr_col, 'p':p_col,
                    'snp':snp_col,'z':z_col,
                    'pos':position_col,'weight':weight,
                    'a1':meta_a, 'a2':a2, 'pvalcol':FAM_P_TITLE,'markercol':FAM_SNP_TITLE}
    print("Index Dict:")
    print fam_indices
    return fam_indices

def read_yank_titles(line):
    title_list = line.strip().split()
    lz_col = title_list.index(YANK_LZ_TITLE)
    chr_col = title_list.index(YANK_CHR_TITLE)
    p_col = title_list.index(YANK_P_TITLE)
    im_col = title_list.index(YANK_IM_TITLE)
    pos_col = title_list.index(YANK_POS_TITLE)
    start_col = title_list.index(YANK_REG_START_TITLE)
    end_col = title_list.index(YANK_REG_END_TITLE)
    id_col = title_list.index(YANK_BAND_ID_TITLE)
    yank_indices = {'lz':lz_col,'chr':chr_col,'start':start_col,
                    'end':end_col,'pos':pos_col,'im':im_col,'p':p_col,'id':id_col}
    #print("Index Dict:")
    #print yank_indices
    return yank_indices

def log_folder(output_folder):
    log_folder = os.path.join(output_folder, 'logs')
    if not os.path.exists(log_folder):
        os.makedirs(log_folder)
    return log_folder

def chr_folder(output_folder, chromosome):
    '''Creates an output subfolder based on the output folder and the chromosome number.
        Args:
            output_folder -- filepath of output mother folder
            chromosome -- chromosome number in string form
    '''
    output_subfolder = os.path.join(output_folder, 'chr'+str(chromosome))
    if not os.path.exists(output_subfolder):
        os.makedirs(output_subfolder)
    return output_subfolder

def read_dict(dict_loc):
    new_dict = dict()
    if dict_loc is not None:
        with open(dict_loc, mode = "r") as dict_file:
            for line in dict_file:
                line_list = line.strip().split()
                new_dict[line_list[0]]=line_list[1]
    return new_dict

def retrieve_freq(freq_loc, snp_name):
    maf_col = 4
    a1_col = 2
    a2_col = 3
    snp_col = 1
    with open(freq_loc, mode = "r") as freq_file:
##        print("Now in retrieve_freq")
        for line in freq_file:
            line_split = line.strip().split()
            if line_split[snp_col]==snp_name:
                return (line_split[maf_col],line_split[a1_col],line_split[a2_col])

def retrieve_r2(snpstar, snpi, ld_loc):
    line1 = True
    with open(ld_loc, mode="r") as  ld:
        for line in ld:
            if line1:
                line1=False
            else:
                lsplit = line.strip().split()
                if lsplit[2] == snpstar and lsplit[5]==snpi:
                    return lsplit[6]

def retrieve_im(map_loc, snp_name, repair_loc=None):
    pos_col = 3
    name_col = 1
    chr_col = 0
    #repair_dict = dict()
    #if repair_loc is not None:
    repair_dict = read_dict(repair_loc)
    if snp_name in repair_dict:
        return repair_dict[snp_name]
    elif snp_name.startswith('rs'):
        return snp_name
    else:
        chrnum = snp_name.partition(':')[0]
        pos = snp_name.partition(':')[2]
        ch = chrnum[3:]
        
        with open(map_loc, mode="r") as map_file:
            for line in map_file:
                split = line.strip().split()
                if split[chr_col]==ch and split[pos_col]==pos:
                    return split[name_col]
                
def is_number(s):
    '''checks a string to determine if it is a number

    '''
    try:
        float(s)
        return True
    except ValueError:
        return False

    
def correct_snp(index_dict, line_list):
    snp = line_list[index_dict['snp']]
    pos = line_list[index_dict['pos']]
    ch = line_list[index_dict['chr']]
    if snp.startswith('rs'):
        return snp , snp
    else:
        lz_snp = 'chr'+ch+':'+pos
        return snp, lz_snp

def in_interval(position, left_endpoint, right_endpoint):
    '''
    Determine if the SNP position falls between the start and end positions of
    the region.

    Args:
        position -- SNP position
        left_endpoint -- start of region
        right_endpoint -- end of region
    '''
    if left_endpoint <=position and position <= right_endpoint:
        return True
    else:
        return False
    
def read_regional_titles(line):
    '''
    Read first line of Region List and assign indices of columns to dictionary.
    Args:
        line -- first line of region list
    Returns:
        reg_indices - dictionary of: 'start','end','sym' indices
    '''
    title_list = line.strip().split()
    reg_chr_col = title_list.index(REG_CHR_TITLE)
    start_col = title_list.index(START_TITLE)
    end_col = title_list.index(END_TITLE)
    gene_name_col = title_list.index(GENE_TITLE)
    band_col = title_list.index(BAND_TITLE)
    title_col = title_list.index(TITLE_TITLE)
    id_col = title_list.index(ID_TITLE)
    reg_indices = {'chr':reg_chr_col,'start':start_col,'end':end_col,
                   'sym':gene_name_col,'band':band_col,
                   'title':title_col,'ID':id_col}
    return reg_indices


def correct_gene_sym(gene_symbol, chromosome):
    ''' Determine if reference gene symbol is appropriate for dictionary use,
        and change where appropriate.
        Args:
            gene_symbol -- gene_symbol listed in the region list
            chromosome -- chromosome listed in region list
        Returns:
            new_symbol -- symbol to be used in dictionary
    '''
    new_symbol = gene_symbol
    if gene_symbol=='0' or gene_symbol =='No' or gene_symbol == 'no' or gene_symbol == 'NA':
        new_symbol = 'nogene-'+ chromosome
    elif gene_symbol=='multiple' or gene_symbol=='Multiple':
        new_symbol = 'multiple-'+chromosome
    return new_symbol


def create_region_list(region_loc):
    region_list = list()
    #RegInfo = namedtuple('RegInfo', 'chro,start,end,sym,band,title,ID')
    standardizer = 1
##    if position_form == 'mb':
##        standardizer = 1e6
##    elif position_form == 'kb':
##        standardizer = 1e3
    with open(region_loc, mode = "r") as region_file:
        line1 = True
        for line in region_file:
            if line.rstrip() == "":
                continue
            if line1:
                reg_indices = read_regional_titles(line)
                line1 = False
            else:
                line_split = line.strip().split()
                chromosome = str(line_split[reg_indices['chr']])
                start_str = line_split[reg_indices['start']]
                end_str = line_split[reg_indices['end']]
                band_name = line_split[reg_indices['band']]
                title_name = line_split[reg_indices['title']]
                id_name = line_split[reg_indices['ID']]
##                start_float = float(line_split[reg_indices['start']])
##                start_int = int(start_float*standardizer)
##                end_int = int(float(line_split[reg_indices['end']])*standardizer)
                #RegInfo = namedtuple('RegInfo', 'chro,start,end,sym,band,title')
                gene_sym = correct_gene_sym(line_split[reg_indices['sym']],
                                            chromosome)
                region_tup = RegInfo(chro=chromosome,start=start_str,
                                     end=end_str,sym=gene_sym,band=band_name,
                                     title=title_name,ID=id_name)
                #print region_tup
                region_list.append(region_tup)
                print region_tup
    return region_list

def create_bp_region_list(region_loc):
    mb_region_list = create_region_list(region_loc)
    region_list = list()
    for mb_region in mb_region_list:
        region = RegInfo(chro = mb_region.chro, start = int(float(mb_region.start) * 1e6),
                         end = int(float(mb_region.end) *1e6),
                         band=mb_region.band, sym=mb_region.sym,
                         ID=mb_region.ID,title=mb_region.title)
        region_list.append(region)
    return region_list


def assoc_tuple(assoc_line, index_dict):
    '''creates named tuple of the form (chro=chromosome
                                        'p'=p-value,
                                        pos=position
                                        'snp'=SNP,
                                        'or'=odds ratio,
                                        't'=t-statistic,
                                        'lo'=lower confidence interval edge,
                                        'hi' = upper ci edge,
                                        'a1' = A1 allele)

    Keyword arguments:
    line_list -- line split into its list of elements

    Returns - namedtuple

    '''
    line_list = assoc_line.strip().split()
    p_store = line_list[index_dict['p']]
    if is_number(line_list[index_dict['p']] ):
        p_store = float(line_list[index_dict['p']])
    LineInfo = namedtuple('LineInfo', 'p,snp,OR,t,lo,hi,a1,chro,pos')
    log_tuple = LineInfo(p=p_store,
                         snp=line_list[index_dict['snp']],
                         OR=tuple_helper(line_list,index_dict['or']),
                         t=tuple_helper(line_list,index_dict['t']),
                         lo=tuple_helper(line_list,index_dict['lo']),
                         hi=tuple_helper(line_list,index_dict['hi']),
                         a1=tuple_helper(line_list,index_dict['a1']),
                         chro=line_list[index_dict['chr']],
                         pos=line_list[index_dict['pos']])
    return log_tuple

def tuple_helper(line_list, index):
    '''Deals with situations in which there is no index, because there
        is no column containing the specific information.
    '''
    item = ''
    if index is not None:
        item = line_list[int(index)]
    sys.stdout.flush()
    return item


def read_assoc_titles(title_line, c_interval, plink_test='logistic'):
    '''Creates a dictionary holding index locations for relevant information
using the first line of the assocation file

    Keyword arguments:
    title_line --   the first line of data file, containing column titles
    c_interval --   confidence interval, decimal form
    plink_test --   plink test used to produce association file
    Returns - tuple containing the indices of the P and SNP column titles

    '''
    title_list = title_line.strip().split()
    #p_index = first_line_list.index('P')
    #snp_index = first_line_list.index('SNP')
    a1_index = title_list.index('A1')
    #chr_index = first_line_list.index('CHR')
    pos_index = title_list.index('BP')
    or_index = None
    stat_index = None
    ci_hi_index = None
    ci_lo_index = None
    nmiss_index = None
    se_index = None
    #THE NEXT 2 TRYS ARE ME BEING LAZY, THEY SHOULD BE CHECKED AND ADDED TO IFS
    try:
        nmiss_index = title_list.index('NMISS')
    except ValueError:
        nmiss_index = None
    try:
        se_index = title_list.index('SE')
    except ValueError:
        se_index = None
    if plink_test in ('logistic', 'fisher','assoc', 'linear'):
        try:
            or_index = title_list.index('OR')
        except ValueError:
            or_index = None
    if plink_test in ('linear', 'assoc'):
        or_index = title_list.index('BETA')
    if plink_test in ('logistic','linear'):
        stat_index = title_list.index('STAT')
    if c_interval is not None:
        ci = str(int(float(c_interval)*100))
        ci_lo_index = title_list.index('L'+ci)
        ci_hi_index = title_list.index('U'+ci)
    index_dict = {'p':title_list.index('P'),
                  'chr':title_list.index('CHR'),
                  'pos':pos_index,
                  'snp':title_list.index('SNP'),
                  'or':or_index,
                  't':stat_index,
                  'hi':ci_hi_index,
                  'lo':ci_lo_index,
                  'a1':a1_index,
                  'nmiss':nmiss_index,
                  'se':se_index,
                  'pvalcol':'P',
                  'markercol':'SNP'}
    return index_dict

def identify_map_indices():
    index_dict = {'chr':0,'pos':3,'snp':1}
    return index_dict

def read_meta_titles(line):
    title_list = line.strip().split()
    meta_chr_col = title_list.index(META_CHR_TITLE)
    p_col = title_list.index(P_TITLE)
    snp_col = title_list.index(SNP_TITLE)
    position_col = title_list.index(POSITION_TITLE)
    z_col = title_list.index(Z_TITLE)
    weight = title_list.index(WEIGHT_TITLE)
    meta_a = title_list.index(A1_TITLE)
    a2 = title_list.index(A2_TITLE)
    direction_col = title_list.index(DIRECTION_TITLE)
    meta_indices = {'chr':meta_chr_col, 'p':p_col,
                    'snp':snp_col,'z':z_col,
                    'pos':position_col,'weight':weight,
                    'a1':meta_a, 'a2':a2,'direction':direction_col,
                    'pvalcol':P_TITLE,
                    'markercol':SNP_TITLE}
    return meta_indices

def create_annot_dict(build, purpose):
    annot_dict_loc = fix_it.locate_annot_dict(build)
    annot_dict = fix_it.build_annot_dict(purpose,annot_dict_loc)
    return annot_dict

def give_table_annotation(table_loc,new_table_loc,index_dict, annot_dict):
    new = open(new_table_loc, mode = "w")
    line1 = True
    with open(table_loc, mode="r") as table:
        for line in table:
            line_list = line.strip().split()
            if line1:
                line_list.append('annotation')
                line1 = False
            else:
                snp = line_list[index_dict['snp']]
                annot = fix_it.get_lz_annot(annot_dict,snp)
                line_list.append(annot)
            new.write('\t'.join(line_list)+'\n')
    
