#! /usr/bin/python2.7
#! ./
'''
Created on Dec 14, 2011

@author Jessica Bonnie
'''
from collections import namedtuple

REG_CHR_TITLE = 'gene_chr'
START_TITLE = 'region_start'
END_TITLE = 'region_end'
GENE_TITLE = 'gene_symbol'
BAND_TITLE = 'Chr_band'




def retrieve_freq(freq_loc, snp_name):
    maf_col = 4
    a1_col = 2
    a2_col = 3
    snp_col = 1
    with open(freq_loc, mode = "r") as freq_file:
        for line in freq_file:
            line_split = line.strip().split()
            if line_split[snp_col]==snp_name:
                return (line_split[maf_col],line_split[a1_col],line_split[a2_col])

def retrieve_im(map_loc, snp_name):
    pos_col = 3
    name_col = 1
    chr_col = 0
    if snp_name.startswith('rs'):
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
    reg_indices = {'chr':reg_chr_col,'start':start_col,'end':end_col,
                   'sym':gene_name_col,'band':band_col}
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
    region_list = []
    RegInfo = namedtuple('RegInfo', 'chro,start,end,sym')
    standardizer = 1
##    if position_form == 'mb':
##        standardizer = 1e6
##    elif position_form == 'kb':
##        standardizer = 1e3
    with open(region_loc, mode = "r") as region_file:
        line1 = True
        for line in region_file:
            if line1:
                reg_indices = read_regional_titles(line)
                line1 = False
            else:
                line_split = line.strip().split()
                chromosome = str(line_split[reg_indices['chr']])
                start_str = line_split[reg_indices['start']]
                end_str = line_split[reg_indices['end']]
                band_name = line_split[reg_indices['band']]
##                start_float = float(line_split[reg_indices['start']])
##                start_int = int(start_float*standardizer)
##                end_int = int(float(line_split[reg_indices['end']])*standardizer)
                RegInfo = namedtuple('RegInfo', 'chro,start,end,sym,band')
                gene_sym = correct_gene_sym(line_split[reg_indices['sym']],
                                            chromosome)
                region_tup = RegInfo(chro=chromosome,start=start_str,
                                     end=end_str,sym=gene_sym,band=band_name)
                #print region_tup
                region_list.append(region_tup)
    return region_list

##def create_region_dict(region_loc, position_form):
##    '''
##    Read file containing regional information and create a dictionary holding
##    that information.
##    Args:
##        region_loc -- filepath of file containing region information
##        position_form -- string ('kb' or 'mb')describing units used in region file
##    Returns:
##        gene_region_dict -- dictionary of region information
##    
##    '''
##    #create a dictionary to hold the regional information
##    gene_region_dict = {}
##    region_list = []
##    RegInfo = namedtuple('RegInfo', 'chro,start,end,sym')
##    
##    standardizer = 1
##    if position_form == 'mb':
##        standardizer = 1e6
##    elif position_form == 'kb':
##        standardizer = 1e3
##    
##    with open(region_loc, mode='r') as region_file:
##        line1 = True
##        for line in region_file:
##            if line1:
##                reg_indices = read_regional_titles(line)
##                line1 = False
##            else:
##                
##                line_list = line.strip().split()
##                gene_sym = line_list[reg_indices['sym']]
##                chr_num = line_list[reg_indices['chr']]
##                new_sym= correct_gene_sym(gene_sym,chr_num)
##                start_fo = int(float(line_list[reg_indices['start']])*standardizer)
##                end_fo = int(float(line_list[reg_indices['end']])*standardizer)
##                reg_tuple = RegInfo(chro=chr_num, start=start_fo, end=end_fo,sym=new_sym)
##                gene_region_dict.update({new_sym
##                                         :[int(chr_num),
##                                           start_fo,
##                                           end_fo]})
##    return gene_region_dict

def assoc_tuple(line_list, index_dict):
    '''creates named tuple of the form ('p'=p-value,
                                        'snp'=SNP,
                                        'or'=odds ratio,
                                        't'=t-statistic,
                                        'lo'=lower confidence interval edge,
                                        'hi' = upper ci edge,
                                        'a1' = A1 allele)

    Keyword arguments:
    line_list -- line split into its list of elements

    Returns - tuple of the form (p-value, SNP, odds ratio, t-statistic)

    '''
    p_store = line_list[index_dict['p']]
    if is_number(line_list[index_dict['p']] ):
        p_store = float(line_list[index_dict['p']])
    LineInfo = namedtuple('LineInfo', 'p,snp,OR,t,lo,hi,a1')
    log_tuple = LineInfo(p=p_store,
                         snp=line_list[index_dict['snp']],
                         OR=tuple_helper(line_list,index_dict['or']),
                         t=tuple_helper(line_list,index_dict['t']),
                         lo=tuple_helper(line_list,index_dict['lo']),
                         hi=tuple_helper(line_list,index_dict['hi']),
                         a1=tuple_helper(line_list,index_dict['a1']))
    return log_tuple

def identify_assoc_cols(first_line_list, c_interval, plink_test='logistic'):
    '''assigns locations of P and SNP columns to global variables

    Keyword arguments:
    first_line_list -- list created from the first line of data file,
        containing column titles
    Returns - tuple containing the indices of the P and SNP column titles

    '''
    
    p_index = first_line_list.index('P')
    snp_index = first_line_list.index('SNP')
    a1_index = first_line_list.index('A1')
    or_index = None
    stat_index = None
    ci_hi_index = None
    ci_lo_index = None
    if plink_test in ('logistic', 'fisher','assoc', 'linear'):
        or_index = first_line_list.index('OR')
    if plink_test in ('linear', 'assoc'):
        or_index = first_line_list.index('BETA')
    if plink_test in ('logistic','linear'):
        stat_index = first_line_list.index('STAT')
    if c_interval is not None:
        ci = str(int(float(c_interval)*100))
        ci_lo_index = first_line_list.index('L'+ci)
        ci_hi_index = first_line_list.index('U'+ci)
    index_dict = {'p':p_index,
                  'snp':snp_index,
                  'or':or_index,
                  't':stat_index,
                  'hi':ci_lo_index,
                  'lo':ci_hi_index,
                  'a1':a1_index}
    return index_dict

