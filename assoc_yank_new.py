#! /ust/bin/python2.7
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


global region_loc
##global position_form
global assoc_loc
global out_file


##
##REG_CHR_TITLE = 'gene_chr'
##START_TITLE = 'region_start'
##END_TITLE = 'region_end'
##GENE_TITLE = 'gene_symbol'
##
##
##ASSOC_CHR_TITLE = 'CHR'
##SNP_TITLE = 'SNP'
##P_TITLE = 'P'
##POSITION_TITLE = 'BP'
##STAT = 'STAT'
##OR = 'OR'
##A1 = 'A1'

C_INTERVAL = .95

#Defines path location of allele frequency file
FREQ_LOC = '/home/jkb4y/work/data/UK_12212011/UK_control.freq'

#Defines path location of original, un map_adapted map file
MAP_LOC = '/home/jkb4y/work/data/UK_12212011/copy_of_UK.bim'

ASSOC_LOC = '/home/jkb4y/work/Projects/AA/data/AA.assoc.logistic'
REGION_LOC ='/home/jkb4y/work/Projects/AA/data/Region_Lists/T1D_regions.txt'
OUT_FILE = '/home/jkb4y/work/results/UK_12212011/yank_results/T1D_assocLeastPs.txt'


##def in_interval(position, left_endpoint, right_endpoint):
##    '''
##    Determine if the SNP position falls between the start and end positions of
##    the region.
##
##    Args:
##        position -- SNP position
##        left_endpoint -- start of region
##        right_endpoint -- end of region
##    '''
##    if left_endpoint <=position and position <= right_endpoint:
##        return True
##    else:
##        return False

def usage():
    print('''
USAGE: assoc_yank.py [FLAG] OBJECT
      FLAG                  DESCRIPTION                                 CURRENT DEFAULT
    -o, --outfile           pathname of  output file                    {0}
    -a, --assoc             pathname of association file                {1}
    -f, --freq              pathname of frequency file                  {2}
    -m, --map               path of unadapted map file                  {3}
    -r, --region-list       pathname of region list file                {4}
    -h, --help              display this usage string



'''.format(OUT_FILE, ASSOC_LOC, FREQ_LOC, MAP_LOC, REGION_LOC))
    
##def read_regional_titles(line):
##    '''
##    Read first line of Region List and assign indices of columns to globals.
##    Args:
##        line -- first line of region list
##    '''
##    title_list = line.strip().split()
##    reg_chr_col = title_list.index(REG_CHR_TITLE)
##    start_col = title_list.index(START_TITLE)
##    end_col = title_list.index(END_TITLE)
##    gene_name_col = title_list.index(GENE_TITLE)
##    reg_indices = {'chr':reg_chr_col,'start':start_col,'end':end_col,'sym':gene_name_col}
##    return reg_indices

##def correct_gene_sym(gene_symbol, chromosome):
##    ''' Determine if reference gene symbol is appropriate for dictionary use,
##        and change where appropriate.
##        Args:
##            gene_symbol -- gene_symbol listed in the region list
##            chromosome -- chromosome listed in region list
##        Returns:
##            new_symbol -- symbol to be used in dictionary
##    
##    '''
##    new_symbol = gene_symbol
##    if gene_symbol=='0' or gene_symbol =='No' or gene_symbol == 'no' or gene_symbol == 'NA':
##        new_symbol = 'nogene-'+ chromosome
##    elif gene_symbol=='multiple' or gene_symbol=='Multiple':
##        new_symbol = 'multiple-'+chromosome
##    return new_symbol

##    
##
##def read_assoc_titles(line):
##    title_list = line.split()
##    chr_col = title_list.index(ASSOC_CHR_TITLE)
##    p_col = title_list.index(P_TITLE)
##    snp_col = title_list.index(SNP_TITLE)
##    pos_col = title_list.index(POSITION_TITLE)
##    or_col = title_list.index(OR)
##    stat_col = title_list.index(STAT)
##    a1_col = title_list.index(A1)
##    if c_interval is not None:
##        ci = str(int(float(c_interval)*100))
##        ci_lo_index = title_list.index('L'+ci)
##        ci_hi_index = title_list.index('U'+ci)
##    assoc_indices = {'chr':chr_col, 'p':p_col,'snp':snp_col,'pos':pos_col,
##                     'stat':stat_col,'OR':or_col,'a1':a1_col}
##    return assoc_indices
    
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
##                reg_indices = pc_toolbox.read_regional_titles(line)
##                line1 = False
##            else:
##                line_list = line.strip().split()
##                gene_sym = line_list[reg_indices['sym']]
##                chr_num = line_list[reg_indices['chr']]
##                line_list[reg_indices['sym']]=correct_gene_sym(gene_sym,chr_num)
##                gene_region_dict.update({line_list[reg_indices['sym']]
##                                         :[int(chr_num),
##                                           int(float(line_list[reg_indices['start']])*standardizer),
##                                           int(float(line_list[reg_indices['end']])*standardizer)]})
##    return gene_region_dict

def read_meta(assoc_loc, gene):
    line1 = True
    standardizer = 1e6
##    if position_form == 'mb':
##        standardizer = 1e6
##    elif position_form == 'kb':
##        standardizer = 1e3
    'or':or_index,
                  't':stat_index,
                  'hi':ci_lo_index,
                  'lo':ci_hi_index,
                  'a1':a1_index,
    old_p = 1
    old_out={}
    old_t = 0
    reg_chr = gene.chro
    reg_start = int(float(gene.start) *standardizer)
    reg_end = int(float(gene.end) * standardizer)
    sym = gene.sym
    band = gene.band
    assoc_counter = 0
    multi_counter = 1
    zero_list = []
    with open(assoc_loc, mode="r") as assoc:
        for line in meta:
            if line1:
                  assoc_indices = pc_toolbox.read_assoc_titles(line)
                  assoc_counter = 1
                  line1 = False
                  old_out = {'chr':'','band':'','start':'',
                             'end':'','sym':'','snp':'',
                             'p':'1','pos':'', 'or':'0',
                             'hi':'0', 'lo':'0','a1':'--',
                             't':'0','note':'--','checkme':'0'}
            else:
                line_info = read_assoc(line,.95, 'logistic')
##
##
##
##
##NOTE NOTE NOTE THIS IS WHERE EDITING STOPPED!!!! RESUME CHANGES HERE!!!!!!!!!!!!!
##
##
##
##                
                line_list = line.strip().split()
                snp_pos = int(line_list[assoc_indices['pos']])
                if int(reg_chr) == int(line_list[assoc_indices['chr']]
                                  ) and pc_toolbox.in_interval(snp_pos, reg_start, reg_end):
                    
                    old_p = float(old_out['p'])
                    cur_p = float(line_list[assoc_indices['p']])
                    cur_snp = line_list[assoc_indices['snp']]
                    cur_t = float(line_list[assoc_indices['t']])
                    cur_hi = line_list[assoc_indices['hi']]
                    cur_lo = line_list[assoc_indices['lo']]
                    cur_a1 = line_list[assoc_indices['a1']]
                    cur_or = line_list[assoc_indices['or']]
##                    if not cur_snp.startswith('rs'):
##                        lz_snp = 'chr'+str(reg_chr)+':'+ str(snp_pos)
##                    else:
##                        lz_snp = cur_snp
##                    if cur_w > weight_min or (reg_chr in ['23','24'] and cur_w > sex_weight):
##                        cur_out = {'chr':gene.chro,'band':gene.band,'start':gene.start,
##                                   'end':gene.end,'sym':gene.sym,'snp':cur_snp,
##                                   'lz':lz_snp,'p':str(cur_p),'pos':str(snp_pos),
##                                   '|z|':str(cur_abs_z),'z':cur_z,'note':'--','checkme':'0',
##                                   'weight':str(cur_w), 'meta_a1':cur_a1,'meta_a2':cur_a2}
##                    else:
                    cur_out = {'chr':gene.chro,'band':gene.band,'start':gene.start,
                               'end':gene.end,'sym':gene.sym,'snp':cur_snp,
                               'p':str(cur_p),'pos':str(snp_pos),'t':cur_t,
                               'note':'--','or':cur_or,
                               'checkme':'1','hi':cur_hi, 'lo':cur_lo,'a1':cur_a1}
                    if cur_p == old_p:
                        multi_counter = multi_counter + 1
                        ZTup = namedtuple('ZTup','snp,p,t,OR')
                        cur_tup = ZTup(snp=cur_out['snp'],OR=cur_out['or'],
                                       p=cur_out['p'],t=cur_out['t'])
                        old_tup = ZTup(snp=old_out['snp'],OR=cur_out['or'],
                                       p=old_out['p'],t=old_out['t'])
                        zero_list.append(cur_tup)
                        if old_tup not in zero_list:
                            zero_list.append(old_tup)
##                        print('NOTE: Both {0} and {1} have p-values of {2}.'.format(old_out['snp'],cur_snp, cur_p))
##                        print('{0} has a z-score of: {1}'.format(old_out['snp'],old_out['z']))
##                        print('{0} has a z-score of: {1}'.format(cur_out['snp'],cur_out['z']))
##                        if cur_out['|z|'] > old_out['|z|']:
##                            old_p = cur_p
##                            old_abs_z = cur_abs_z
##                            old_out.update(cur_out)
##                        print('Retaining {0} as the significant SNP based on z-score value.'.format(old_out['snp']))
##                        old_out['note']='{0}SNPs(p={1})'.format(multi_counter, old_p)
                    elif cur_p < old_p:
                        multi_counter = 1
                        old_p = cur_p
                        old_out.update(cur_out)
                        zero_list = []
                meta_counter = meta_counter + 1
    sorted_list = sorted(zero_list, key=lambda member: member[3], reverse=True)
    print old_out
    return old_out, sorted_list




def read_assoc(assoc_loc, gene_name,gene_region_dict):
    line1 = True
    
    old_p = 1
    old_out=[]
    reg_chr, reg_start, reg_end = gene_region_dict[gene_name]
    assoc_counter = 0
    multi_counter = 1
    sig_list = []
    blank_dict = {'chr':'--','start':'--','end':'--','sym':'--','snp':'--',
                  'lz':'--','pval':'--','stat':'--','OR':'NA','aa1':'NA',
                  'pos':'--','notes':'--','maf':'NA','ma1':'NA','ma2':'NA'}
    order_list = ['chr','start','end','sym','snp','lz','pval','stat','OR',
                  'aa1','maf','ma1','ma2','notes']
    with open(assoc_loc, mode="r") as assoc:
        for line in assoc:
            line_list = line.strip().split()
            if line1:
                assoc_i = pc_toolbox.read_assoc_titles(line, c_interval)
##                assoc_indices = read_assoc_titles(line)
                assoc_counter = 1
                line1 = False
            else:
                snp_pos = int(line_list[assoc_i['pos']])
                if not line_list[assoc_i['p']] == 'NA':
                    if reg_chr == int(line_list[assoc_i['chr']]
                                      ) and pc_toolbox.in_interval(snp_pos,reg_start, reg_end):
                        cur_p = float(line_list[assoc_i['p']])
                        cur_snp , lz_snp = pc_toolbox.correct_snp(assoc_i,line_list)
##                        cur_snp = line_list[assoc_indices['snp']]
##                        lz_snp = '--'
##                        if not cur_snp.startswith('rs'):
##                            lz_snp = 'chr'+str(reg_chr)+':'+ str(snp_pos)
##                        else:
##                            lz_snp = cur_snp
                        cur_out = [str(reg_chr), str(reg_start), str(reg_end),gene_name,
                                   cur_snp,lz_snp, str(cur_p), str(snp_pos),'--']
                        if cur_p == old_p:
                            multi_counter = multi_counter + 1
                            cur_tup = (cur_snp, lz_snp,str(cur_p))
                            old_tup = (old_out[4], old_out[5],str(old_p))
                            sig_list.append(cur_tup)
                            if old_tup not in sig_list:
                                sig_list.append(old_tup)
                            print('NOTE: Both {0} and {1} have p-values of {2}.'.format(old_out[4],cur_snp, cur_p))
                            print('There are no Z-scores with which to decide between them, so I arbitrarily retain {0}.'.format(old_out[4]))
                            old_out[8]='{0}SNPs(p={1})'.format(multi_counter, old_p)
                            print(old_out[8])
                        elif cur_p < old_p:
                            multi_counter = 1
                            old_p = cur_p
                            old_out = cur_out
                            sig_list = []
                    assoc_counter = assoc_counter + 1
    
    return old_out, sig_list


def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global assoc_loc, region_loc, out_file, freq_loc,c_interval
    
    position_form = POSITION_FORM
    assoc_loc = ASSOC_LOC
    region_loc = REGION_LOC
    out_file = OUT_FILE
    freq_loc = FREQ_LOC
    map_loc = MAP_LOC
    c_interval = C_INTERVAL
    
    try: 
        opts, args = getopt.getopt(argv, "ha:r:o:m:",
                                   ["help","assoc=","region-list=",
                                    "out=","map=","ci="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
##        elif opt in ("--bp-form"):
##            position_form = arg
        elif opt in ("-a","--assoc"):
            assoc_loc = arg
        elif opt in ("-f","--freq"):
            freq_loc = arg
        elif opt in ("-m","--map"):
            map_loc = arg
        elif opt in ("-r","--region-list"):
            region_loc = arg
        elif opt in ("-o","--out"):
            out_file = arg
        elif opt in ("--ci"):
            c_interval = arg
        
def find_leastP_SNPs(gene_region_list, assoc_loc,out_file, freq_loc):
    global out_file
    title_list = ['chr','chr_band','region_start(Mb)',
                  'region_end(Mb)','gene_symbol',
                  'snp_name','snp_position','locuszoom_snp',
                  'p-value','notes','OR','t-stat','MAF',
                  'MAF_A1','MAF_A2']
    sig_title_list = ['snp_name','locuszoom_snp','p-value','MAF','MAF_A1','MAF_A2']
    title_line = '\t'.join(title_list)
    sig_line = '\t'.join(sig_title_list)
    with open(out_file, mode="w") as out_text:
        out_text.write(title_line +'\n')
    with open(out_file, mode="a") as out_text:
        for gene in gene_list:
            gene_info, sig_list = read_assoc(assoc_loc, gene,gene_region_dict)
            if not sig_list == []:
                sig_start, ext = os.path.splitext(out_file)
                sig_path = sig_start+'_'+'Chr'+gene_info[0] + '_'+ gene_info[3]+'.txt'
                sigs = open(sig_path,mode='w')
                sigs.write(sig_line +'\n')
                for sig in sig_list:
                    sigs.write('\t'.join(sig)+'\n')
                sigs.close()
            gene_info[1]=str(float(gene_info[1])*1e-6)
            gene_info[2]=str(float(gene_info[2])*1e-6)
            gene_line = '\t'.join(gene_info[:9])
            print(gene_info[:9])
            out_text.write(gene_line + '\n')
            create_condition_list(gene_info)
            
def create_condition_list(gene_info):
    global out_file
    if gene_info[8]== '--':
        
        head, tail = os.path.split(out_file)
        subfolder_path = os.path.join(head, 'assoc_condition_lists')
        if not os.path.exists(subfolder_path):
            os.makedirs(subfolder_path)
        snp_name = gene_info[4]
        chr_num = gene_info[0]
        position = gene_info[7]
        gene_sym = gene_info[3]
        mb_range = gene_info[1] + '-'+gene_info[2]
        filestart = os.path.join(subfolder_path,'Chr'+chr_num + '_')
        filename = filestart + gene_sym + '.txt'
        clist = open(filename, mode='w')
        clist.write(gene_info[5])
        clist.close()
    
    
def main():
    global out_file, region_loc,assoc_loc, freq_loc
    cl_arguments(sys.argv[1:])
    gene_region_list = pc_toolbox.create_region_list(region_loc)
    #gene_region_dict = create_region_dict(region_loc,position_form)
    #key_list = gene_region_dict.keys()
    #gene_list = sorted(key_list,key=lambda gene:gene_region_dict[gene][0])
    #find_leastP_SNPs(gene_list, gene_region_dict, assoc_loc)
    find_leastP_SNPs(gene_region_list, assoc_loc, out_file, freq_loc)
                

if __name__ == '__main__':
    main()
    
