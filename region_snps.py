#! /usr/bin/python2.6
#! ./
'''
Created Oct 28, 2011

@author: Jessica Bonnie
'''
import os
import re
import sys
import pc_toolbox
import fix_it

BUILD = 'hg18'
REGION_LOC = '/home/jkb4y/work/data/Region_Lists/hg18/T1D_region_PTPN2.txt'
META_LOC = '/home/jkb4y/work/data/2012Feb1/eurmeta/eurmeta_03062012.txt'
WEIGHT_MIN = 21000
OUTFOLDER = '/home/jkb4y/work/results/eurmeta/'





def read_meta(meta_loc, region, outfolder,build, annot_dict=None):
    line1 = True
    standardizer = 1e6
    reg_chr = region.chro
    reg_start = int(float(region.start) *standardizer)
    reg_end = int(float(region.end) * standardizer)
    sym = region.sym
    band = region.band
    outfile = os.path.join(outfolder,band+'.tbl')
    out = open(outfile, mode = "w")
    title_list = ['imchip_name','rs_name','chr','{0}_pos'.format(build),'p-value','z-score','weight','a1','a2']
    out.write('\t'.join(title_list)+'\n')
    with open(meta_loc, mode="r") as meta:
        for line in meta:
            if line1:
                meta_indices = pc_toolbox.read_meta_titles(line)
                line1 = False
            else:
                line_list = line.strip().split()
                snp_pos = int(line_list[meta_indices['pos']])
                weight = line_list[meta_indices['weight']]
                if int(reg_chr) == int(line_list[meta_indices['chr']]
                                       ) and pc_toolbox.in_interval(snp_pos, reg_start, reg_end):
                    if float(weight)>WEIGHT_MIN:
                        p_val = line_list[meta_indices['p']]
                        lz_snp = line_list[meta_indices['snp']]
                        z_score = line_list[meta_indices['z']]
                        a1 = line_list[meta_indices['a1']]
                        a2 = line_list[meta_indices['a2']]
                        snp_rs = annot_dict[lz_snp].rs
                        snp_im = annot_dict[lz_snp].name
                        out.write('\t'.join([snp_im,snp_rs,reg_chr,str(snp_pos),p_val,z_score,weight,a1,a2])+'\n')
    out.close()
                        
                    




def main(argv):
    region_loc = REGION_LOC
    build = BUILD
    meta_loc = META_LOC
    outfolder = OUTFOLDER
    region_list = pc_toolbox.create_region_list(region_loc)
    annot_dict_loc = fix_it.locate_annot_dict(build)
    annot_dict = fix_it.build_annot_dict('LOG',annot_dict_loc)
    fixed_meta_loc = fix_it.locate_fixed_meta(meta_loc, build)
    for region in region_list:
        read_meta(fixed_meta_loc, region, outfolder, build, annot_dict)




if __name__ == '__main__':
    main(sys.argv[1:])
