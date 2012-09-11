#! /usr/bin/python2.7
#! ./
'''
Created on Nov 2, 2011

@author: Jessica Bonnie
'''
import os
import subprocess
import sys
import getopt
from collections import namedtuple
import pc_toolbox


global table_loc
##global fixed_table_loc
global outfolder
global region_loc
global assoc
global population
global fix

##ASSOC_CHRCOL = 0
##META_CHRCOL = 7
##ASSOC_POSCOL = 2
##POSCOL = 8
##ASSOC_MARKERCOL = 1
##MARKERCOL = 0
##ASSOC_PVALCOL = 8
##PVALCOL =  5
#BAD_TABLE_LOC = '/home/jkb4y/work/data/UK_12212011/meta.tbl'
TABLE_LOC = '/home/jkb4y/work/data/2011Nov28/aaedit.assoc.logistic'
OUTFOLDER = '/home/jkb4y/work/results/2011Nov28/AssocGraphs-YRI/'
OUT_FLAG = ''
CHROMOSOME = None
REGION_LOC = None
DATA_SET = None
MAX_BP = 5e7
POPULATION = 'CEU'
ASSOC = False

RegInfo = namedtuple('RegInfo', 'chro,start,end,sym,band')



def new_table_loc(table_loc):
    global assoc
    base, ext = os.path.splitext(table_loc)
    if assoc is not None:
        base, dext = os.path.splitext(base)
        ext = dext+ext
    new_base = base + '_lz'
    lz_table_loc = new_base + ext
    return lz_table_loc


def edit_table(table_loc, fix):
    line1=True
    ztup_list = []
    lz_table_loc = new_table_loc(table_loc)
    with open(table_loc, mode = 'r') as table:
        if fix:
            lz_table = open(lz_table_loc, mode='w')
        #with open(lz_table_loc, mode='w')as lz_table:
        for line in table:
            line_split = line.strip().split()
            if line1:
                if assoc:
                    index_dict = pc_toolbox.read_assoc_titles(line, c_interval=None,plink_test='logistic')
                else:
                    index_dict = pc_toolbox.read_meta_titles(line)
                line1=False
            else:               
##                if line1:
##                    table_title_list = line_split
##                    parse_col()
##                    new_table.write('\t'.join(line_split)+'\n')
##                    line1 = False
##                else:
##                    cur_chr = line_split[meta_chrcol]
                cur_chr = line_split[index_dict['chr']]
                cur_pos = int(line_split[index_dict['pos']])
                if not line_split[index_dict['snp']].startswith('rs'):
                    line_split[index_dict['snp']]='chr'+cur_chr+':'+str(cur_pos)
                if not line_split[index_dict['p']]=='NA':
                    if float(line_split[index_dict['p']])==0:
                        ZTup = namedtuple('ZTup','lz,chro,pos')
                        z = ZTup(lz=line_split[index_dict['snp']],
                                 chro=line_split[index_dict['chr']],
                                 pos=line_split[index_dict['pos']])
                        ztup_list.append(z)
            if fix:
                lz_table.write('\t'.join(line_split)+'\n')
    if fix:
        lz_table.close()
        table_loc = lz_table_loc
    return index_dict, ztup_list, table_loc



def read_table(table_loc, chromosome, index_dict):

    #global high_pos, low_pos
    low_pos = 5000000000
    high_pos = 0
    line1 = True
    with open(table_loc, mode='r') as table:
        for line in table:
            line_split = line.strip().split()
            if line1:
                if assoc:
                    index_dict = pc_toolbox.read_assoc_titles(line, c_interval=None,plink_test='logistic')
                else:
                    index_dict = pc_toolbox.read_meta_titles(line)
                line1=False
            else:
                cur_chr = line_split[index_dict['chr']]
                cur_pos = int(line_split[index_dict['pos']])
                if chromosome==cur_chr and cur_pos < low_pos:
                    low_pos = cur_pos
                if chromosome==cur_chr and cur_pos > high_pos:
                    high_pos = cur_pos
    return (low_pos,high_pos)
                    
def macro_lz_title(chromosome, counter, total):
    global assoc
    if assoc:
        graph_title = 'Chromosome {0}: Macroview ({1} of {2})'.format(chromosome, counter, total)
    else:
        graph_title = 'Chromosome {0}: Meta-Analysis Macroview ({1} of {2})'.format(chromosome,counter,total)
    return graph_title

def locuszoom(region,graph_title, index_dict):
    
    global outfolder, population
    global table_loc, out_flag
    markertitle = str(index_dict['snp'])
    pvaltitle = str(index_dict['p'])
    ld_cache = os.path.join(outfolder, 'chr{0}_ld_cache.db'.format(region.chro))
    lz_cl_args = ['locuszoom','--metal', table_loc,
                  '--chr',region.chro, '--start',str(region.start),
                  '--end',str(region.end),'--pvalcol',index_dict['pvalcol'],
                  '--markercol',index_dict['markercol'],
                  '--delim','whitespace', '--cache',ld_cache,
                  '--prefix', region.band+out_flag,
                  '--plotonly','--verbose','--no-date',
                  '--snpset','None',
                  'showRecomb=False',
                  'title='+graph_title, 'ylab=-log(p-value)',
                  '--pop',population]
    p = subprocess.Popen(lz_cl_args, bufsize = 0, executable=None,stdin=None,
                         stdout=None,stderr=None, preexec_fn=None,close_fds=False,
                         shell=False,cwd=None,env=None, universal_newlines=False,
                         startupinfo=None, creationflags=0)
    p.wait()

def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global table_loc,chromosome, fix
    global outfolder, out_flag, region_loc
    global assoc, data_set, max_bp, population

    chromosome = CHROMOSOME
    table_loc = TABLE_LOC
    outfolder = OUTFOLDER
    out_flag = OUT_FLAG
    region_loc = REGION_LOC
    data_set = DATA_SET
    assoc = ASSOC
    max_bp=MAX_BP
    population = POPULATION
    fix = False
    
    
    try: 
        opts, args = getopt.getopt(argv, "ho:f:c:t:r:",
                                   ["help","outfolder=","flag=",
                                    "chromosome=","table=",
                                    "region-list=","pop=",
                                    "data=","max-bp=","assoc","fix",
                                    "region="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-f","--flag"):
            out_flag = arg
        elif opt in ("-o","--outfolder"):
            outfolder = arg
        elif opt in ("-c","--chromosome"):
            chromosome = arg
##        elif opt in ("--meta-chrcol"):
##            meta_col_dict['chr'] = arg
##        elif opt in ("--poscol"):
##            meta_col_dict['pos'] = arg
##        elif opt in ("--markercol"):
##            meta_col_dict['snp'] = arg
        elif opt in ("--max-bp"):
            max_bp= float(arg)
##            meta_col_dict['p'] = arg
        elif opt in ("-t","--table"):
            table_loc = arg
        elif opt in ("-r","--region-list"):
            region_loc = arg
        elif opt in ("--pop"):
            population = arg
            print("Population argument has been read.")
        elif opt in ("--assoc"):
            #meta_col_dict = assoc_col_dict
            assoc = True
        elif opt in ("--fix"):
            #meta_col_dict = assoc_col_dict
            fix = True
        elif opt in ("--data"):
            data_set = arg
        elif opt in ("--region"):
            reg_info = eval(arg)

def usage():
    print('''Usage goes here....        
USAGE: lz_meta.py [FLAG] OBJECT
      FLAG                  DESCRIPTION                                     CURRENT DEFAULT
    -o, --outfolder         path of folder where results are to be written  {0}
    -c, --chromosome        chromosome number for which to draw macroview   {1}
    -t, --table             path to data table                              {2}
    -r, --region-list       path to region list                             {3}
    --max-bp                maximum number of bp contained in a macroview plot {4}
    --assoc
    --pop                   LD population                                   {5}
    --fix                   used to indicate that SNP names in table must be changed for lz
    -f, --flag              flag to add to output file names                {6}
    -h, --help              display this usage string
    '''.format(OUTFOLDER,CHROMOSOME,TABLE_LOC, REGION_LOC, MAX_BP, POPULATION,
               OUT_FLAG))

def determine_range(range_list, start_pos, end_pos):
    global max_bp
    difference = end_pos - start_pos
    if difference < max_bp:
        range_tuple = (start_pos, end_pos)
        #print(range_tuple)
        range_list.append(range_tuple)
    else:
        halfway = int(round((start_pos + (difference/2)),0))
        determine_range(range_list,start_pos,halfway)
        determine_range(range_list,halfway,end_pos)
    return range_list

def macro_ranges(table_loc, chromosome, index_dict):
    (low_pos,high_pos) = read_table(table_loc, chromosome, index_dict)
    print('The base pair range of the snps is from {0} to {1}.'.format(low_pos,high_pos))
    range_list = []
    region_list = []
    unsorted_ranges = determine_range(range_list,low_pos,high_pos)
    ranges = sorted(unsorted_ranges, key=lambda r_tuple:r_tuple[0])
    return ranges
    #return range_list

def macro_regions(chromosome, ranges):
    regions = []
    for couplet in ranges:
        start = int(couplet[0]- int(2e3))
        end = int(couplet[1]+ int(2e3))
        start_mb = round(float(start)*1e-6,1)
        end_mb = round(float(end)*1e-6,1)
        flag = 'Chr'+chromosome+'_'+str(start_mb)+'-'+str(end_mb)
        RegInfo = namedtuple('RegInfo', 'chro,start,end,sym,band')
        reg_info = RegInfo(band=flag,chro=chromosome,start=couplet[0],end=couplet[1],sym=flag)
##        quadlet = (flag, chromosome, couplet[0],couplet[1])
        regions.append(reg_info)
    return regions

def remove_zerops(region_list,zero_list):
    for zero in zero_list:
##        print('Current ZTup being considered:')
##        print zero
        for region in region_list:
##            print('Current region being considered:')
##            print region
            if float(zero.chro) == float(region.chro) and pc_toolbox.in_interval(float(zero.pos),
                                                                                 float(region.start),
                                                                                 float(region.end)):
                print('''
***********************************************************
    
        NOTE: REMOVING Region {0} from list due to:
            {1} at chr{2}:{3} with p=0.
        
***********************************************************
'''.format(region.band,zero.lz,zero.chro,zero.pos))
                sys.stdout.flush()
                region_list.remove(region)
                
    return region_list



def which_one(region_loc, table_loc,chromosome, index_dict, zero_list):
    if region_loc is not None:
##        determine_table_titles(table_loc)
##        region_list = create_region_list(region_loc, region_col_dict)
        mb_region_list = pc_toolbox.create_region_list(region_loc)
        region_list = []
        for mb_region in mb_region_list:
            region = RegInfo(chro = mb_region.chro, start = int(float(mb_region.start) * 1e6),
                             end = int(float(mb_region.end) *1e6), band=mb_region.band, sym=mb_region.sym)
            region_list.append(region)
        print region_list
    else:
        range_list = macro_ranges(table_loc, chromosome, index_dict)
##        print(range_list)
        region_list = macro_regions(chromosome, range_list)
    fix_list = remove_zerops(region_list, zero_list)
    print('Here are the regions which will be drawn:')
    print fix_list
    return fix_list

def correct_gene_symbol(gene_symbol, chromosome):
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

def draw_things(table_loc, region_loc, chromosome, population, index_dict,z_list):
    pop_tag = ''
    if population is not None:
        pop_tag = '\n HapMap {0} LD Population'.format(population)
    
    counter = 1
    region_list = which_one(region_loc, table_loc,chromosome, index_dict, z_list)
    total = len(region_list)
    for region in region_list:
##        chromosome = quadlet[1]
##        flag = quadlet[0]
        if region_loc is None:
            graph_title = macro_lz_title(region.chro, counter, total)+pop_tag
        else:
            #graph_title = 'Region {0}: Meta-Analysis'.format(flag)
            graph_title = 'Region {0}: Analysis'.format(region.band) + pop_tag
        locuszoom(region, graph_title, index_dict)
##        locuszoom(region.chro, str(region.start), str(region.end),region.band,graph_title)
        counter = counter + 1

def main():
    global outfolder, table_loc, chromosome, region_loc
    global data_set,assoc, population, fix
    cl_arguments(sys.argv[1:])
    print("Population is: {0}.".format(population))
    #lz_table = new_table_loc(table_loc)
    #if bad_table_loc is not None:
    index_dict, zero_list, table_loc = edit_table(table_loc, fix)
    os.chdir(outfolder)
    if chromosome is None and region_loc is None:
        for i in range(1,26):
            draw_things(table_loc, region_loc, str(i), population, index_dict, zero_list)
    else:
        draw_things(table_loc, region_loc, chromosome, population, index_dict, zero_list)

if __name__=='__main__':
    main()
