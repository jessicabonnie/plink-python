
import pc_toolbox
import fix_it
import sys
import subprocess
import os

def make_ref_dict(yank_loc):
    ref_dict = dict()
    with open(yank_loc, mode="r") as yank:
        line1 = True
        
        for line in yank:
            line_list = line.strip().split()
            if line1:
                yank_indices = pc_toolbox.read_yank_titles(line)
                line1 = False
            else:
                ref_dict[line_list[yank_indices["id"]]] = (line_list[yank_indices["lz"]],line_list[yank_indices["im"]])
    return ref_dict


def locuszoom(lz_cmd, outfolder):
    print(lz_cmd)
    os.chdir(outfolder)
    p = subprocess.Popen(lz_cmd, bufsize = 0, executable=None,stdin=None,
                         stdout=None,stderr=None, preexec_fn=None,close_fds=False,
                         shell=False,cwd=None,env=None, universal_newlines=False,
                         startupinfo=None, creationflags=0)
    sys.stdout.flush()
    p.wait()
    sys.stdout.flush()

def main(argv):
    global altref, altannot
    altref = False
    altannot = False
    yank_loc = "/home/jkb4y/cphgdesk_share/Projects/IMCHIP/Intersect_SNP_list/eurmeta_10172012/RegionYank/eurmeta_yank.tbl"
    #region_loc = '/home/jkb4y/work/data/Region_Lists/hg19/T1D_regions_hg19_09202012.txt'
    #region_loc = '/home/jkb4y/work/data/Region_Lists/hg19/T1D_regions_hg19_09202012_chr19_20_21_22_23.txt'
    region_loc = '/home/jkb4y/work/data/Region_Lists/hg19/achilleas_all_09202012.txt'
    #region_loc = '/home/jkb4y/work/data/Region_Lists/hg19/achilleas_chr21_22_23_09202012.txt'
    #region_loc = '/home/jkb4y/work/data/Region_Lists/hg19/MHC_zoom.txt'
    #assoc_loc = '/home/jkb4y/ubs/work/results/July2012/ModelComparison/AA_cov_12_nosex.assoc.logistic'
    assoc_loc = '/home/jkb4y/ubs/work/results/July2012/ModelComparison/AA_cov_12_nosex_SIlz.assoc.logistic'

    #outfolder = '/home/jkb4y/ubs/work/results/July2012/RegionGraphs/RefTopUK'
    #outfolder = '/home/jkb4y/cphgdesk_share/Projects/IMCHIP/July2012/RegionGraphs/MHC'
    #outfolder = '/home/jkb4y/ubs/work/results/July2012/RegionGraphs/RefTopUK_expand_0.2'
    #outfolder = '/home/jkb4y/ubs/work/results/July2012/RegionGraphs/SI_Annot_expansion_0.1'
    outfolder = '/home/jkb4y/ubs/work/results/July2012/RegionGraphs/Achilleas_regions'
    outfolder = '/home/jkb4y/cphgdesk_share/Projects/IMCHIP/July2012/RegionGraphs/Achilleas_regions_expansion_0.1'
    
    #assoc_loc = fix_it.locate_fixed_table(assoc_loc, build="hg19",table_type = "assoc")
    regions = pc_toolbox.create_bp_region_list(region_loc)
    expansion = int(.1e6)
    ref_dict = make_ref_dict(yank_loc)
    print ref_dict
    for region in regions:
        prefix = region.ID + '_AA'
        cache_folder = os.path.join(outfolder, 'ld_caches')
        if not os.path.exists(cache_folder):
            os.makedirs(cache_folder)
        ld_cache = os.path.join(cache_folder, 'chr{0}_ld_cache.db'.format(region.chro))
        lz_args = ["locuszoom", "--metal", assoc_loc,
                   '--chr',region.chro, '--start',str(region.start - expansion),
                   '--end',str(region.end + expansion),'--pvalcol','P','--markercol','EXTREME_LZ',
                   "--pop","AFR",'--build','hg19','--source','1000G_Nov2010',
                   '--prefix', prefix,
                   #'annotCol=annotation','annotPch=24,25,22,23,21',
                   #'annotOrder=splice,nonsyn,coding,utr,no-annotation',
                   '--plotonly','--verbose','--no-date','--delim','whitespace','ymax=5',
                   '--cache',ld_cache]
        if altannot:
            lz_args.extend([#'prelude=./altannot_prelude.R',
                            'postlude=/home/jkb4y/h4t1/programs/plink_python/aa_postlude.R',
                            'refsnpTextColor=transparent',
                            'annotCol=SI_annotation',
                            'annotPch=8,21','annotOrder=SI,none'])
        else:
            lz_args.extend(['annotCol=annotation','annotPch=24,25,22,23,21',
                            'annotOrder=splice,nonsyn,coding,utr,no-annotation'])
##        if refsnp in ['rs75793288'] and not altref:
##            lz_args.extend(['--no-ld'])
        #if region.title in ['20p13'] and not altref:
         #   lz_args.extend(['--no-ld'])
        if region.ID in ['1p21.2_y','1q23.3_b_y','12p13.33_y'] and not altref:
            lz_args.extend(['--no-ld'])
        if int(region.chro) > 22:
            lz_args.extend(['--no-ld'])
        if region.title not in ['6p21.32'] and altref:
            refsnp = ref_dict[region.ID][0]
            if refsnp in ['rs75793288','rs3216815','rs12416116']:
                lz_args.extend(['--no-ld'])
            title = "{0}: AA population, refsnp {1}".format(region.ID, refsnp)
            lz_args.extend(["--refsnp",refsnp,'title={0}'.format(title)])
        else:
            lz_args.extend(['title={0}'.format(region.ID)])
            if region.ID in ['20p13']:
                lz_args.extend(['--no-ld'])
            if region.ID in ['6p21.32']:
                lz_args.extend(['rfrows=7'])
        locuszoom(lz_args, outfolder)        

if __name__=='__main__':
    main(sys.argv[1:])
