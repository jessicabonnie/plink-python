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
from operator import attrgetter
import pc_toolbox
import achilleas_yank


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
TABLE_LOC = None
OUTFOLDER = None
OUT_FLAG = ''
CHROMOSOME = None
REGION_LOC = None
DATA_SET = None
MAX_BP = 5e7
POPULATION = 'CEU'
ASSOC = False
BUILD = 'hg19'
CHR_POSTLUDE = '/home/jkb4y/h4t1/programs/plink_python/chr_postlude.R'
POSTLUDE = '/home/jkb4y/h4t1/programs/plink_python/postlude.R'
EQTL_POSTLUDE = '/home/jkb4y/h4t1/programs/plink_python/eqtl_postlude.R'
AA_POSTLUDE = '/home/jkb4y/h4t1/programs/plink_python/aa_postlude.R'
RegInfo = pc_toolbox.RegInfo
MANUSCRIPT_NO_LD_LIST = ['4q27','12q13.3','10q23.31','Xp22.2','Xq28']
AA_NO_LD_LIST = ['20p13','Xp22.2','Xq28']
NO_LD_LIST = []
B_NO_LD_LIST = ['1q32.1','2q12.1_z']
CD4_NO_LD_LIST = ['11p15.5','12q13.3']# '10p11.22_z'
CD8_NO_LD_LIST = ['12q13.3']#'12q13.3'
MONO_NO_LD_LIST = ['18p11.21','12q13.3']
NK_NO_LD_LIST = ['12q13.3']#'12q24.12',

#def edit_table(table_loc, fix):
def edit_table(table_loc, table_type):
    global perm
    line1=True
    ztup_list = list()
##    lz_table_loc = new_table_loc(table_loc)
    with open(table_loc, mode = 'r') as table:
        for line in table:
            line_split = line.strip().split()
            if line1:
                if table_type == 'assoc':
                    index_dict = pc_toolbox.read_assoc_titles(line, c_interval=None,plink_test='logistic')
                elif table_type == 'eqtl':
                    index_dict = achilleas_yank.read_eqtl_titles(line, perm)
                    print index_dict
                else:
                    index_dict = pc_toolbox.read_meta_titles(line)
                    print index_dict
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
##                if not line_split[index_dict['snp']].startswith('rs'):
##                    line_split[index_dict['snp']]='chr'+cur_chr+':'+str(cur_pos)
                if not line_split[index_dict['p']]=='NA':
                    if float(line_split[index_dict['p']])==0:
                        ZTup = namedtuple('ZTup','lz,chro,pos,abs_z')
                        if assoc:
                            z_or_t=abs(float(line_split[index_dict['t']]))
                        else:
                            z_or_t=abs(float(line_split[index_dict['z']]))
                        z = ZTup(lz=line_split[index_dict['snp']],
                                 chro=line_split[index_dict['chr']],
                                 pos=line_split[index_dict['pos']],
                                 abs_z=z_or_t)
                        ztup_list.append(z)
    ztup_sort = sorted(ztup_list, key=lambda z: z.abs_z, reverse=True)
    return index_dict, ztup_sort, table_loc



def read_table(table_loc, chromosome, index_dict):

    #global high_pos, low_pos
    low_pos = 5000000000
    high_pos = 0
    line1 = True
    with open(table_loc, mode='r') as table:
        for line in table:
            line_split = line.strip().split()
            if line1:
##                if assoc:
##                    index_dict = pc_toolbox.read_assoc_titles(line, c_interval=None,plink_test='logistic')
##                elif eqtl:
##                    index_dict = pc_toolbox.read_assoc_titles(line, c_interval=None,plink_test='linear')
##                else:
##                    index_dict = pc_toolbox.read_meta_titles(line)
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
    global table_type, assoc
    #if assoc:
    if table_type in ['assoc','eqtl']:
        graph_title = 'Chromosome {0}: Macroview ({1} of {2})'.format(chromosome, counter, total)
    else:
        graph_title = 'Chromosome {0}: Meta-Analysis Macroview ({1} of {2})'.format(chromosome,counter,total)
    return graph_title


def locuszoom(table_loc, region,graph_title, index_dict, hitspec_loc=None):
    
    global outfolder, population, cell_type
    global out_flag, build, chrom, table_type

    prelude = '/home/jkb4y/h4t1/programs/achilleas_plots/gene_prelude.R'
    no_ld_list = []
    postlude = POSTLUDE
    perm_flag = ''
    if check_perm(table_loc):
        perm_flag = '_perm'
    if chrom:
        postlude = CHR_POSTLUDE
    elif table_type == 'eqtl':
        postlude = EQTL_POSTLUDE
        if out_flag == '':
            out_flag = '_'+cell_type + perm_flag
        if cell_type == 'B':
            no_ld_list = B_NO_LD_LIST
        elif cell_type == 'CD4':
            no_ld_list = CD4_NO_LD_LIST
        elif cell_type == 'CD8':
            no_ld_list = CD8_NO_LD_LIST
        elif cell_type == 'MONO':
            no_ld_list = MONO_NO_LD_LIST
        elif cell_type == 'NK':
            no_ld_list = NK_NO_LD_LIST
    elif population == 'AFR':
        postlude = AA_POSTLUDE
        no_ld_list = AA_NO_LD_LIST
    print('''
                    NOW ENTERING LOCUSZOOM
''')
    b_arg = list()
    if build == 'hg19':
        b_arg.extend(['--build','hg19','--source','1000G_Nov2010','--pop',population])
    if build == 'hg18':
        b_arg.extend(['--build','hg18','--pop',population])
        if region.chro not in ['23','24','25','26']:
            b_arg.extend(['--source','1000G_June2010'])
        #else:
        #b_arg.extend(['--build','hg18','--source','1000G_June2010','--pop',population])#, 'showRecomb=False'])
    sys.stdout.flush()
    markertitle = str(index_dict['snp'])
    pvaltitle = str(index_dict['p'])
    path_dir = os.path.abspath(os.path.join(os.path.dirname(outfolder),".."))
    #path_dir = os.path.dirname(outfolder)
    cache_folder = os.path.join(path_dir, 'ld_caches')
    if not os.path.exists(cache_folder):
        os.makedirs(cache_folder)
    ld_cache = os.path.join(cache_folder, 'chr{0}_ld_cache.db'.format(region.chro))
    lz_cl_args = ['locuszoom','--metal', table_loc,
                  '--chr',region.chro, '--start',str(region.start),
                  '--end',str(region.end),'--pvalcol',index_dict['pvalcol'],
                  '--markercol',index_dict['markercol'],
                  '--delim','whitespace', '--cache',ld_cache,
                  '--prefix', region.band+out_flag,
                  '--plotonly',
                  '--verbose','--no-date',
                  'geneFontSize=.6',
                  '--snpset','None']
    lz_cl_args.extend(b_arg)
    
    if chrom:
        lz_cl_args.extend(['refsnpTextColor=transparent',
                           'rfrows=0',
                           'showRecomb=FALSE',
                           '--no-ld',
                           'ymax=17',
                           'postlude='+postlude,
                           'title=Chromosome {0}'.format(region.chro).replace('23','X')])
    elif table_type == 'eqtl':
        eqtl_title = '{0}_{1}'.format(region.title, cell_type)
        if check_perm(table_loc):
            eqtl_title = eqtl_title + ' (perm)'
        if cell_type == 'CD4':
            if region.ID in ['16q23.1','9p24.2']:
                lz_cl_args.extend(['legend=left'])
        if cell_type == 'B':
            if region.ID in ['3p21.31']:
                lz_cl_args.extend(['legend=left'])
            if region.ID in ['17q12']:
                lz_cl_args.extend(['prelude='+prelude])
            if region.ID in ['6p21.32']:
                lz_cl_args.extend(['legend=right'])
        if cell_type == 'MONO':
            if region.ID in ['7p15.2']:
                lz_cl_args.extend(['legend=right'])
        if cell_type == 'NK':
            if region.ID in ['16p13.13_a']:
                lz_cl_args.extend(['legend=left'])
                
            
        lz_cl_args.extend(['refsnpTextColor=transparent',
                           'ymax=4',
                           'postlude='+postlude,
                           #'prelude='+prelude,
                           'title={0}'.format(eqtl_title)])
        ##options used for manuscript regional graphs:
    else:
        if region.ID == '6p21.32':
            lz_cl_args.extend(['rfrows=10'])
        else:
            lz_cl_args.extend(['rfrows=4'])
        if region.ID == '13q32.3':
            lz_cl_args.extend(['legend=left'])
        lz_cl_args.extend(['ymax=4',#'refsnpTextColor=transparent',
                           'postlude='+postlude,
                           'title={0}'.format(region.title)])
        if region.ID in NO_LD_LIST:
            lz_cl_args.extend(['--no-ld'])
            

        ##options used in hg19 graphs of ANY type:
    if build == 'hg19':
        lz_cl_args.extend(['annotCol=annotation','annotPch=24,25,22,23,21',
                           'annotOrder=splice,nonsyn,coding,utr,no-annotation'])
    if region.ID in no_ld_list:
        lz_cl_args.extend(['--no-ld'])
        
        ##options used for regional graphs of regions which fail using LD:
##                  '--no-ld',
        ##options used for non-manuscript regional graphs:
                  #'title='+graph_title]
    
    print lz_cl_args
##    lz_cl_args.extend(b_arg)
    sys.stdout.flush()
    if hitspec_loc is not None:
        graph_title = graph_title.replace('<','l.t').replace('(','-').replace(')','')
        print graph_title
        title_lines = graph_title.splitlines()
        print title_lines
        title = ' '.join(title_lines)
        print title
        hitspec_list = ['NA',str(region.chro),str(region.start),
                        str(region.end),'NA','yes',
                        "title='"+title+"' snpset=NULL weightCol='Weight' theme='publication'"]
        with open(hitspec_loc,mode="a") as hitty:
            hitty.write('\t'.join(hitspec_list)+'\n')
        return

    p = subprocess.Popen(lz_cl_args, bufsize = 0, executable=None,stdin=None,
                         stdout=None,stderr=None, preexec_fn=None,close_fds=False,
                         shell=False,cwd=None,env=None, universal_newlines=False,
                         startupinfo=None, creationflags=0)
    sys.stdout.flush()
    p.wait()
    sys.stdout.flush()




def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global table_loc,chromosome, build, table_type
    global outfolder, out_flag, region_loc, chrom, cell_type
    global assoc, data_set, max_bp, population, hitspec

    chromosome = CHROMOSOME
    table_loc = TABLE_LOC
    outfolder = OUTFOLDER
    out_flag = OUT_FLAG
    region_loc = REGION_LOC
    data_set = DATA_SET
    assoc = ASSOC
    max_bp=MAX_BP
    population = POPULATION
    build = BUILD
    hitspec = False
    chrom = False
    table_type = 'meta'
    cell_type = ''
    
    
    try: 
        opts, args = getopt.getopt(argv, "ho:f:c:t:r:",
                                   ["help","outfolder=","flag=",
                                    "chromosome=","table=",
                                    "region-list=","pop=",
                                    "data=","max-bp=","assoc","build=",
                                    "region=","hitspec","chrom","eqtl="])
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
            table_type = 'assoc'
            assoc = True
        elif opt in ("--eqtl"):
            #meta_col_dict = assoc_col_dict
            table_type = 'eqtl'
            cell_type = arg
        elif opt in ("--build"):
            #meta_col_dict = assoc_col_dict
            build = arg
##        elif opt in ("--data"):
##            data_set = arg
##        elif opt in ("--region"):
##            reg_info = eval(arg)
        elif opt in ("--hitspec"):
            hitspec = True
        elif opt in ("--chrom"):
            chrom = True

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
    --build                 genome build: hg18 or hg19                      {6}
    -f, --flag              flag to add to output file names                {7}
    -h, --help              display this usage string
    '''.format(OUTFOLDER,CHROMOSOME,TABLE_LOC, REGION_LOC, MAX_BP, POPULATION,BUILD,
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
    print('''
*********************************************************************
    The base pair range of the snps is from {0} to {1}.
*********************************************************************
'''.format(low_pos,high_pos))
    sys.stdout.flush()
    range_list = list()
    #region_list = list()
    unsorted_ranges = determine_range(range_list,low_pos,high_pos)
    ranges = sorted(unsorted_ranges, key=lambda r_tuple:r_tuple[0])
    return ranges
    #return range_list

def macro_regions(chromosome, ranges):
    regions = list()
    for couplet in ranges:
        start = int(couplet[0]- int(2e3))
        end = int(couplet[1]+ int(2e3))
        start_mb = round(float(start)*1e-6,2)
        end_mb = round(float(end)*1e-6,2)
        flag = 'Chr'+chromosome+'_'+str(start_mb)+'-'+str(end_mb)
        reg_info = RegInfo(band=flag,chro=chromosome,start=couplet[0],end=couplet[1],sym=flag,title=flag,ID=flag)
        quadlet = (flag, chromosome, couplet[0],couplet[1])
        regions.append(reg_info)
    return regions

def identify_zerops(region_list,zero_list):
    z_region_list = list()
    for zero in zero_list:
        #print('Current ZTup being considered:')
        #print zero
        for region in region_list:
##            print('Current Region being considered:')
##            print region
            if str(zero.chro) == str(region.chro) and pc_toolbox.in_interval(float(zero.pos),
                                                                                 float(region.start),
                                                                                 float(region.end)):
                if region not in z_region_list:
                    z_region_list.append(region)
                    print('''
***********************************************************
    
        NOTE: REMOVING Region {0} from first list due to:
            {1} at chr{2}:{3} with p=0.
        
***********************************************************
'''.format(region.ID,zero.lz,zero.chro,zero.pos))
                    sys.stdout.flush()
                    continue
                #region_list.remove(region)
    sys.stdout.flush()
    return region_list, z_region_list



def which_one(region_loc, table_loc,chromosome, index_dict, zero_list):
    if region_loc is not None:
##        determine_table_titles(table_loc)
##        region_list = create_region_list(region_loc, region_col_dict)
        region_list = pc_toolbox.create_bp_region_list(region_loc)
##        mb_region_list = pc_toolbox.create_region_list(region_loc)
##        region_list = list()
##        for mb_region in mb_region_list:
##            region = RegInfo(chro = mb_region.chro, start = int(float(mb_region.start) * 1e6),
##                             end = int(float(mb_region.end) *1e6),
##                             band=mb_region.band, sym=mb_region.sym,
##                             ID=mb_region.ID,title=mb_region.title)
##            region_list.append(region)
        print region_list
    else:
        range_list = macro_ranges(table_loc, chromosome, index_dict)
##        print(range_list)
        region_list = macro_regions(chromosome, range_list)
##    base, ext = os.path.splitext(table_loc)
##    noZ_loc = base + '_noZs'+ext
    fix_list, zr_list = identify_zerops(region_list, zero_list)
    print('''
************************************************************************************
Here are the regions which will be drawn:
''')
    for fix in fix_list:
        print fix
    print('''
The following regions contain zero p-values, and will be adapted:''')
    for z in zr_list:
        print z
        
    print('''
************************************************************************************''')
    sys.stdout.flush()
    return fix_list, zr_list

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
    global build, hitspec
    hitspec_loc = None
    if hitspec:
        base, ext = os.path.splitext(table_loc)
        hitspec_loc = base + '_hitspec.txt'
        with open( hitspec_loc,mode="w") as hitty:
            hitty.write('\t'.join(['snp','chr','start','stop',
                                   'flank','run','m2zargs'])+'\n')
    pop_tag = ''
##    if build == 'hg18':
##        pop_tag = '\n 1000G(June2010) {0} LD Population'.format(population)
##    if build == 'hg19':
##        pop_tag = '\n 1000G(Nov2010) {0} LD Population'.format(population)
    
    base, ext = os.path.splitext(table_loc)
    noZ_loc = base + '_noZs'+ext
    noZ_loc = table_loc
    
    counter = 1
    region_list, zr_list = which_one(region_loc, noZ_loc,chromosome, index_dict, z_list)
##    region_list, zr_list = which_one(region_loc, table_loc,chromosome, index_dict, z_list)
    total = len(region_list)
    for region in region_list:
        if build == 'hg18':
            if region.chro in ['23','24','25','26']:
                pop_tag = '\n HAPMAP {0} LD Population'.format(population)
            else:
                pop_tag = '\n 1000G(June2010) {0} LD Population'.format(population)
##        chromosome = quadlet[1]
##        flag = quadlet[0]
        if region_loc is None:
            graph_title = macro_lz_title(region.chro, counter, total)
        else:
            #graph_title = 'Region {0}: Meta-Analysis'.format(flag)
            graph_title = '{0}'.format(region.title)
            
        if region not in zr_list:
            #t_loc = table_loc
            graph_title = graph_title + pop_tag
        else:
            #t_loc = noZ_loc
            p_tag ='\nReference SNP P-Value < 1e-100'
##            if int(region.chro) in [6, 11, 1]:
##                p_tag ='\n(0->1e-100)'
##            else:
##                p_tag ='\n(0->1e-65)'
            graph_title = graph_title + p_tag
            
        sys.stdout.flush()
        locuszoom(noZ_loc, region, graph_title, index_dict, hitspec_loc)
##        locuszoom(t_loc, region, graph_title, index_dict)
        sys.stdout.flush()
##        locuszoom(region.chro, str(region.start), str(region.end),region.band,graph_title)
        counter = counter + 1
##    counter = 1
##    for z in zr_list:
##        print z.chro
##        p_tag = ''
##        if int(z.chro) == 6:
##            p_tag ='\n(0->1e-350)'
##        else:
##            p_tag ='\n(0->1e-150)'
##        if region_loc is None:
##            g_title = macro_lz_title(z.chro, counter, total)
##        else:
##            g_title = 'Region {0}: Analysis'.format(z.band)
##        graph_title = g_title + p_tag
##        sys.stdout.flush()
##        locuszoom(noZ_loc, z, graph_title, index_dict)
##        sys.stdout.flush()
##        counter = counter + 1
def check_perm(table_loc):
    answer = False
    if 'perm' in table_loc:
        answer = True
    return answer

def main(argv):
    print("Arguments to lz_meta.py:")
    print argv
    global outfolder, table_loc, chromosome, region_loc, cell_type
    global data_set,assoc, population, build, hitspec, table_type, perm
    cl_arguments(argv)
    print('''
************************************************************************************
LocusZoom instructed that population is: {0}.
************************************************************************************'''.format(population))
    sys.stdout.flush()
    #lz_table = new_table_loc(table_loc)
    #if bad_table_loc is not None:
    perm = check_perm(table_loc)
    index_dict, zero_list, table_loc = edit_table(table_loc, table_type)
    sys.stdout.flush()
    os.chdir(outfolder)
    if build == 'hg18':
        if population not in ("YRI","CEU","JPT+CHB"):
            print("ERROR: {0} is not one of the supported populations in hg18!".format(population))
            print("Locuszoom will use LD information from CEU population instead!")
            population = 'CEU'
    if build == 'hg19':
        if population not in ("AFR","EUR","ASN"):
            print("ERROR: {0} is not one of the supported populations in hg19!".format(population))
            print("Locuszoom will use LD information from EUR population instead!")
            population = 'EUR'
    if chromosome is None and region_loc is None:
        for i in range(1,22):
            sys.stdout.flush()
            draw_things(table_loc, region_loc, str(i), population, index_dict, zero_list)
            sys.stdout.flush()
    else:
        sys.stdout.flush()
        draw_things(table_loc, region_loc, chromosome, population, index_dict, zero_list)
        sys.stdout.flush()

if __name__=='__main__':
    main(sys.argv[1:])
