import os
import pc_toolbox
import sys

def filter_by_region(region, file_in, file_out, index_dict):
    line1 = True
    with open(file_in, mode="r") as ifiley:
        with open(file_out, mode="w") as ofiley:
            for line in ifiley:
                line_split = line.strip().split()
                if line1:
                    ofiley.write('\t'.join(line_split)+'\n')
                    line1=False
                else:
                    bp = line_split[index_dict['pos']]
                    chro = line_split[index_dict['chr']]
                    if int(region.chro) == int(chro) and pc_toolbox.in_interval(int(bp),region.start,region.end):
                        ofiley.write('\t'.join(line_split)+'\n')
                        
                    



def main(argv):
    out_folder = '/home/jkb4y/cphgdesk_share/Yan/UBASH3A'
    eurmeta_out = os.path.join(out_folder, 'eurmeta_UBASH3A_122112.txt')
    eurmeta_in = '/home/jkb4y/work/data/2012Feb1/eurmeta_Oct2012/eurmeta_lz_hg19_intersect.tbl'
    eurmeta_dict, thing = pc_toolbox.create_index_dict(eurmeta_in, 'meta')
    family_out = os.path.join(out_folder, 'eurgdtscan_UBASH3A_122112.txt')
    family_in = '/home/jkb4y/work/data/2012Feb1/Family_Nov2012/eurgdtscan_lz_hg19_intersect.txt'
    family_dict, thing = pc_toolbox.create_index_dict(family_in, 'family')
    region_loc ='/home/jkb4y/work/data/Region_Lists/hg19/T1D_regions_hg19_09202012.txt'
    region_list = pc_toolbox.create_bp_region_list(region_loc)
    for region in region_list:
        if region.ID == '21q22.3_a':
            print region
            filter_by_region(region, eurmeta_in,eurmeta_out, eurmeta_dict)
            filter_by_region(region, family_in, family_out, family_dict)
            



if __name__ == '__main__':
    main(sys.argv[1:])    
