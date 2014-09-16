#! /usr/bin/python2.7
#! ./
'''
Created on Nov 2, 2011

@author: Jessica Bonnie
'''
import os
import subprocess
import sys
sys.path.insert(0, '/home/jkb4y/h4t1/programs/plink_python/')
import pc_toolbox
#REGION_LOC = '/h2/t74/cphgdesk/share/cphg_OnengutLab/Jessica/Projects/IMCHIP/Region_Lists/T1D_regions_hg19_04152013_11p15.5.txt'
REGION_LOC = '/h2/t74/cphgdesk/share/cphg_OnengutLab/Jessica/Projects/IMCHIP/Region_Lists/T1D_regions_hg19_04152013.txt'
#REGION_LOC = '/home/jkb4y/work/data/Region_Lists/hg18/T1D_regions_chr1.txt'
BFILE = '/h2/t74/cphgdesk/share/cphg_OnengutLab/Jessica/Projects/IMCHIP/Intersect_SNP_list/2012Oct17/CaseControl/Data/intersect_07232013'
LDFOLDER = '/h2/t74/cphgdesk/share/cphg_OnengutLab/Jessica/Projects/IMCHIP/Intersect_SNP_list/2012Oct17/eurmeta/eurmeta_LD_07232013/'
OUTFOLDER = '/h2/t74/cphgdesk/share/cphg_OnengutLab/Jessica/Projects/IMCHIP/Intersect_SNP_list/2012Oct17/CaseControl/CAnalysis_eurmeta_LD_07232013/'
ASSOC_FOLDER = '/h2/t74/cphgdesk/share/cphg_OnengutLab/Jessica/Projects/IMCHIP/Intersect_SNP_list/2012Oct17/CaseControl/CAnalysis_eurmeta_07232013/'

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

def quick_lz(region, cloop, snp, ld, outfolder, assoc_folder, r2, draw_me=True):
    chrband = region.ID
    chromosome = region.chro
    label = '{0}_~{1}_SNPsOut'.format(chrband,cloop)
    assoc_loc = os.path.join(assoc_folder,label+'.assoc.logistic')
    lz_label = label + '_r2_'+r2
    #lz_label = 'EVI5_EXPERIMENT_'+label + '_r2_'+r2
    title = '{0}:{1} SNPs Out'.format(chrband, cloop)
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
		  lz_label,'--plotonly','--verbose','--no-date',
                  '--snpset','None',
                  #'--cache', ld_cache,
                  #'showRecomb=FALSE','rfrows=1',
                  #'title='+title,
                  '--source','1000G_Nov2010',
                  '--build','hg19','--pop','EUR',
                  'rfrows=4',
                  '--ld',ld]
    if cloop == 0:
        lz_cl_args.extend(['--refsnp',snp])
        title = title + '(LD RefSNP: {0})'.format(snp)
##    else:
##        title = title + '\n(List Provided at Start)'
    title = title + ' (r^2 > {0})'.format(r2)
    lz_cl_args.extend(['title='+title])
    if not draw_me:
        return
    p = subprocess.Popen(lz_cl_args, bufsize = 0, executable=None,stdin=None,
                         stdout=None,stderr=None, preexec_fn=None,close_fds=False,
                         shell=False,cwd=None,env=None, universal_newlines=False,
                         startupinfo=None, creationflags=0)
    p.wait()

def plink_ld(region, bfile, ldfolder, assoc_folder, r2, sum_only=False):
    chromosome = region.chro
    chrband = region.ID
    ldfile = chrband + '_r2_{0}'.format(r2)
    #ldfile = 'EVI5_EXPERIMENT_r2_{0}'.format(r2)
    ld_sub = pc_toolbox.chr_folder(ldfolder, chromosome)
    ld = os.path.join(ld_sub,ldfile)
    #snp_list = 'EVI5_EXPERIMENT.txt'
    snp_list = '{0}_~leastP_SNPs.txt'.format(chrband)
    snp_loc = os.path.join(assoc_folder,snp_list)
    if sum_only:
        return ld + '.ld', snp_loc
    plink_args = ['plink','--noweb','--bfile',bfile, '--chr', chromosome,
                 '--from-mb', region.start, '--to-mb',region.end,
                 '--r2','--ld-window-r2',r2,'--ld-window','999999',
                 '--ld-window-kb','99999','--ld-snp-list',snp_loc,
                 '--out', ld,'--filter-controls','--maf', '0.005']
    p = subprocess.Popen(plink_args, bufsize = 0, executable=None,stdin=None,
                         stdout=None,stderr=None, preexec_fn=None,close_fds=False,
                         shell=False,cwd=None,env=None, universal_newlines=False,
                         startupinfo=None, creationflags=0)
    p.wait()
    return ld + '.ld', snp_loc

def ld_for_lz(ld):
    
    base, ext = os.path.splitext(ld)
    lz_ld = base + '_lz'+ext
    fix = open(lz_ld, mode = "w")
    line1 = True
    with open(ld, mode = "r") as ld_info:
        for line in ld_info:
            if line1:
                fix.write('\t'.join(['snp1','snp2','rsquare','dprime'])+'\n')
                line1 = False
            else:
                split = line.strip().split()
                fix.write(split[5]+'\t'+split[2]+'\t'+split[6]+'\tNA\n')
    fix.close()
    return lz_ld

def copy_to_tot(total_loc, ld_loc):
    line1 = True
    with open(total_loc, mode="a") as total:
        with open(ld_loc, mode="r") as ld:
            for line in ld:
                if line1:
                    line1=False
                    continue
                else:
                    lsplit = line.strip().split()
                    total.write('\t'.join(lsplit)+'\n')
                

def main_lz(region, outfolder, assoc_folder, snp_list,ld, r2, draw_me=True):
    listy = read_list(snp_list)
    print ld
    cloop = 0
    for snp in listy:
        quick_lz(region,cloop,snp, ld, outfolder, assoc_folder, r2, draw_me)
        cloop = cloop + 1
        
def main():
    draw_me = False
    sum_only = False
    MHC = False
    regions = pc_toolbox.create_region_list(REGION_LOC)
    bfile = BFILE
    #ldfolder = LDFOLDER
    ldfolder = '/home/jkb4y/ubs/work/results/Intersection_06022014/eurmeta_LD/'
    #outfolder = OUTFOLDER
    outfolder = '/home/jkb4y/ubs/work/results/Intersection_06022014/CAnalysis_LD/'
    assoc_folder = ASSOC_FOLDER
    #assoc_folder = '/home/jkb4y/work/results/Intersection_06022014/CAnalysis_eurmeta'
    #for r2 in ['0','0.2','0.4','0.6','0.8']:
    for r2 in ['0']:
    #for r2 in ['0','0.2']:
    #for r2 in ['0.2','0.4']:
    #for r2 in ['0.8','0.6']:
    #r2 = '0.8'
        tot_loc = os.path.join(ldfolder,'all_regions_r2_{0}.ld'.format(r2))
        with open(tot_loc, mode="w") as total:
            ld_title_list = ['CHR_A','BP_A','SNP_A','CHR_B','BP_B','SNP_B','R2']
            total.write('\t'.join(ld_title_list)+'\n')
        for region in regions:
            chrband = region.ID
            chromosome = region.chro
            subfolder = pc_toolbox.chr_folder(outfolder, chromosome)
            assoc_chr_folder = pc_toolbox.chr_folder(assoc_folder, chromosome)
            if chrband == '6p21.32' and not MHC:
                continue
            ld, snp_loc = plink_ld(region, bfile, ldfolder, assoc_chr_folder, r2, sum_only)
            copy_to_tot(tot_loc, ld)
            
            lz_ld = ld_for_lz(ld)
            main_lz(region, subfolder, assoc_chr_folder, snp_loc,lz_ld, r2, draw_me)
        

if __name__=='__main__':
    main()
