'''
@author Jessica Bonnie
Created Aug 30, 2012
'''

import os
import sys
import pc_toolbox
import subprocess

def locate_prelude(outfolder, gene, cell_type, region_ID):
    prelude_gene_name = cell_type + '_' + region_ID + '_'+gene +'_'+ "prelude.R"
    prelude_gene_loc = os.path.join(outfolder, prelude_gene_name)
    return prelude_gene_loc

def write_prelude(base_prelude_loc, gene, gene_prelude_loc):
    prelude = open(gene_prelude_loc, mode = "w")
    prelude.write("gene = '"+gene+"'\n")
    base = open(base_prelude_loc, mode = "r")
    for line in base:
        prelude.write(line)

def create_prefix(cell_type, region_ID, gene=None, perm=True, cis=True):
    prefix = region_ID + '_' + cell_type
    if perm:
        prefix = prefix + '_perm'
    if gene is not None:
        prefix = prefix + '_' + gene
    if cis:
        prefix = prefix + '_cis'
    return prefix

def create_title(cell_type, region, gene=None, perm=True):
    title = "{0} Cells: {1}".format(cell_type, region)
##    if perm:
##        title = title + ' (perm)'
    if gene is not None:
        title = title + ', {0}'.format(gene)
    return title

def read_list(list_loc):
    listy = list()
    try:
        with open(list_loc, mode="r") as filey:
            for line in filey:
                lstrip = line.strip()
                if len(lstrip) > 0:
                    listy.append(lstrip)
    except IOError:
        pass
    return listy

def locate_postlude(cell_type, postlude_folder):
    postlude_loc = os.path.join(postlude_folder,'eqtl_postlude_{0}.R'.format(cell_type))
    return postlude_loc

def locuszoom(table_loc, region, graph_title, index_dict, outfolder, refsnp=None, prefix=None, prelude=None, postlude=None):
    
    global population
    global out_flag, build, table_type
    print('''
                    NOW ENTERING LOCUSZOOM
''')
    print(refsnp)
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
    print(index_dict)
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
                  '--markercol','EXTREME_LZ',#index_dict['markercol'],
                  '--delim','whitespace', '--cache',ld_cache,
                  '--prefix', prefix,
                  '--plotonly',
                  '--verbose','--no-date',
                  'geneFontSize=.6',
                  '--snpset','None','ymax=4',
                  'title={0}'.format(graph_title)]
    lz_cl_args.extend(b_arg)
    if refsnp in ['chr2:204576923','chr2:204732714','chr1:206898330',
                  'chr12:58218799','chr6:32549808','chr6:32904980',
                  'chr6:31322486','chr6:32974400']:
        lz_cl_args.extend(['--no-ld'])
        #refsnp = None
##    if refsnp == 'chr21:43823910':
##        refsnp = None
    if refsnp is not None:
        lz_cl_args.extend(['--refsnp',refsnp])
    if postlude is not None:
        lz_cl_args.extend(['refsnpTextColor=transparent',
                           'postlude='+postlude])
    if prelude is not None:
        lz_cl_args.extend(['prelude='+prelude])

        ##options used in hg19 graphs of ANY type:
    if build == 'hg19':
        lz_cl_args.extend(['annotCol=annotation','annotPch=24,25,22,23,21',
                           'annotOrder=splice,nonsyn,coding,utr,no-annotation'])
##    if region.ID in no_ld_list:
##        lz_cl_args.extend(['--no-ld'])
    
    print lz_cl_args
##    lz_cl_args.extend(b_arg)
    sys.stdout.flush()
    p = subprocess.Popen(lz_cl_args, bufsize = 0, executable=None,stdin=None,
                         stdout=None,stderr=None, preexec_fn=None,close_fds=False,
                         shell=False,cwd=None,env=None, universal_newlines=False,
                         startupinfo=None, creationflags=0)
    sys.stdout.flush()
    p.wait()
    sys.stdout.flush()


def read_single(file_loc):
    line1 = True
    with open(file_loc)as filey:
        for line in filey:
            if line1:
                return line.strip()
        

def main(argv):
    global table_type, region_loc, cell_type, table_loc, perm, build, population
    build = 'hg19'
    cis = True
    population = 'EUR'
    region_loc = '/home/jkb4y/work/data/Region_Lists/hg19/T1D_regions_hg19_05242012_extra.5MB.txt'
    prelude_base = '/home/jkb4y/h4t1/programs/achilleas_plots/gene_prelude.R'
    cell_type = "CD8"
    postlude_folder= '/home/jkb4y/h4t1/programs/achilleas_plots/postludes'
    tablefold = '/home/jkb4y/ubs/work/data/Achilleas/cis-eQTL/'
    table_loc = os.path.join(tablefold,'{0}_cis_eqtls_permwBP_JB_LZ.txt'.format(cell_type))
    table_type = 'perm'
    index_dict, thing = pc_toolbox.create_index_dict(table_loc, table_type)
    perm = True
    achilleas_cell_results = "/home/jkb4y/ubs/work/results/Achilleas/hg19/{0}".format(cell_type)
    genefolder = os.path.join(achilleas_cell_results,'GeneLists')
    base_outfolder = os.path.join(achilleas_cell_results, 'PerGene')
    if not os.path.exists(base_outfolder):
        os.makedirs(base_outfolder)
    regions = pc_toolbox.create_bp_region_list(region_loc)
    for region in regions:
        ID = region.ID
        band = region.band
        outfolder = os.path.join(base_outfolder, ID)
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
        os.chdir(outfolder)
        region_genefolder = os.path.join(genefolder, ID)
        genelist = os.path.join(region_genefolder, ID+'_genes.list')
        genes = read_list(genelist)
        for gene in genes:
##                outfolder = os.path.join(base_outfolder, region.ID)
##                if not os.path.exists(outfolder):
##                    os.makedirs(outfolder)
            refsnp_loc = os.path.join(region_genefolder, gene + '_ref.txt')
            refsnp = read_single(refsnp_loc)
            gene_prelude = locate_prelude(outfolder, gene, cell_type, ID)
            cell_postlude = None#locate_postlude(cell_type, postlude_folder)
            write_prelude(prelude_base, gene, gene_prelude)
            graph_title = create_title(cell_type, ID, gene)
            prefix = create_prefix(cell_type, band, gene, perm, cis)
##                os.chdir(outfolder)
            locuszoom(table_loc, region,graph_title, index_dict, outfolder,
                      refsnp, prefix, gene_prelude, cell_postlude)
        os.chdir(base_outfolder)
        try:
            os.rmdir(outfolder)
        except OSError:
            pass



if __name__=='__main__':
    main(sys.argv[1:])
