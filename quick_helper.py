#! /usr/bin/python2.7
#! ./
'''
Created on Nov 2, 2011

@author: Jessica Bonnie
'''
import os
import subprocess


import pc_toolbox
#chr3:46374178
#chr2:162832883

#2	2q24.2	162.6	163.1	IFIH1
REGION_LOC = '/home/jkb4y/work/data/Region_Lists/hg18/T1D_regions_chr1.txt'
BFILE = '/home/jkb4y/work/data/2012Feb1/hg18/imm'
#CHRBAND = '6p21.32'
#ASSOC_LOC = '/home/jkb4y/work/results/2012Feb1/hg18/CAnalysis_test/{0}.assoc.logistic'.format(LABEL)
#NEW_LOC = '/home/jkb4y/work/plink_python/results/UK/chr6/Chr6:MHC_NoZeros.assoc.logistic'
OUT_FOLDER = '/home/jkb4y/work/results/2012Feb1/hg18/CAnalysis_test/'

CHROMOSOME = '6'
MB_START = 24.7
MB_END = 34
#SNP_LIST = '/home/jkb4y/work/results/2012Feb1/hg18/CAnalysis_test/{0}_~leastP_SNPs.txt'.format(CHRBAND)
#LD_FILE = '/home/jkb4y/work/results/2012Feb1/hg18/ld_caches/{0}_fix.ld'.format(CHRBAND)
#1 SNPs Conditioned Out \n(Condition List Provided at Start)'
#0 SNPs Conditioned Out \n({0} Provided as LD Reference SNP'.format(SNP)
#\n(Condition List Provided at Start)'
#METAL_LOC = '/home/jkb4y/work/plink_python/data/meta_lz.tbl'

def read_list(list_loc):
    listy = list()
    with open(list_loc) as list_file:
        for thing in list_file:
            if len(thing.strip().split())> 0:
                listy.append(thing.strip())
    return listy

def fix_assoc(assoc_loc, new_loc):
    line1 = True
    with open(assoc_loc, mode="r") as assoc:
        with open(new_loc, mode = "w") as new:
            for line in assoc:
                line_split = line.strip().split()
                if line1:
                    new.write('\t'.join(line_split)+'\n')
                    line1 = False
                else:
                    if not line_split[8]=='NA':
                        if not float(line_split[8])==0:
                            new.write('\t'.join(line_split)+'\n')

def quick_lz(region, cloop, snp, ld, outfolder):
    chrband = region.band
    chromosome = region.chro
    label = '{0}_~{1}_SNPsOut'.format(chrband,cloop)
    assoc_loc = os.path.join(outfolder,label+'.assoc.logistic')
    title = '{0}: {1} SNPs Conditioned Out'.format(chrband, cloop)
    #title = TITLE
    bp_start = str(int(float(region.start) * 1e6))
    bp_end = str(int(float(region.end) * 1e6))
    os.chdir(outfolder)
    #cwd = os.getcwd()
    #path_dir = os.path.abspath(os.path.join(os.path.dirname(cwd),".."))
    #cache_folder = os.path.join(path_dir, 'ld_caches')
    #ld_cache = os.path.join(cache_folder, 'chr{0}_ld_cache.db'.format(chromosome))
    
    lz_cl_args = ['locuszoom','--metal', str(assoc_loc),
                  '--chr',str(chromosome),'--start', str(bp_start),
                  '--end', str(bp_end), '--pvalcol','P',
                  '--markercol','SNP','--delim','whitespace', '--prefix',
		  label,'--plotonly','--verbose','--no-date',
                  '--snpset','None',
                  #'--cache', ld_cache,
                  #'showRecomb=FALSE','rfrows=1',
                  #'title='+title,
                  'ylab=-log(p-value)',
                  '--source','1000G_June2010',
                  '--build','hg18','--pop','CEU',
                  '--ld',ld]
    if cloop == 0:
        lz_cl_args.extend(['--refsnp',snp])
        title = title + '\n({0} Provided as LD Reference SNP'.format(snp)
    else:
        title = title + '\n(Condition List Provided at Start)'
    lz_cl_args.extend(['title='+title])
##    lz_cl_args = ['locuszoom','--metal', metal_loc,
##                  '--chr',chromosome,'--start', bp_start,
##                  '--end', bp_end, '--pvalcol','P',
##                  '--markercol','SNP','--delim','whitespace', '--prefix',
##		  label,'--plotonly','--verbose','--no-date',
##                  '--snpset','None','showRecomb=FALSE','title='+title]
    p = subprocess.Popen(lz_cl_args, bufsize = 0, executable=None,stdin=None,
                         stdout=None,stderr=None, preexec_fn=None,close_fds=False,
                         shell=False,cwd=None,env=None, universal_newlines=False,
                         startupinfo=None, creationflags=0)
    p.wait()

def plink_ld(chromosome, start, end, bfile, ld_snp_list, out):
    plink_args = ['plink','--noweb','--bfile',bfile, '--chr', chromosome,
                 '--from-mb', start, '--to-mb',end,
                 '--r2','--ld-window-r2','0','--ld-window','999999',
                 '--ld-window-kb','99999','--ld-snp-list',ld_snp_list,
                 '--out', out]
    p = subprocess.Popen(plink_args, bufsize = 0, executable=None,stdin=None,
                         stdout=None,stderr=None, preexec_fn=None,close_fds=False,
                         shell=False,cwd=None,env=None, universal_newlines=False,
                         startupinfo=None, creationflags=0)
    p.wait()

def main_lz(region, outfolder, snp_list,ld):
    listy = read_list(snp_list)
    print ld
    cloop = 0
    for snp in listy:
        quick_lz(region,cloop,snp, ld, outfolder)
        cloop = cloop + 1

def main():
    regions = pc_toolbox.create_region_list(REGION_LOC)
    ldfolder = '/home/jkb4y/work/results/2012Feb1/hg18/LD/'
    outfolder = '/home/jkb4y/work/results/2012Feb1/hg18/CAnalysis_Test/'
    for region in regions:
        chrband = region.band
        chromosome = region.chro
        subfolder = pc_toolbox.chr_folder(outfolder, chromosome)
        ldfile = chrband + '_r2_0.8_lz.ld'
        ld_sub = pc_toolbox.chr_folder(ldfolder, chromosome)
        ld = os.path.join(ld_sub,ldfile)
        snp_list = '{0}_~leastP_SNPs.txt'.format(chrband)
        snp_loc = os.path.join(ld_sub,snp_list)
        main_lz(region, subfolder, snp_loc,ld)
if __name__=='__main__':
    main()
