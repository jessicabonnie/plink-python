#! /usr/bin/python2.6
#! ./
'''
Created Feb 29, 2012

@author: Jessica Bonnie
'''
import os
import getopt
import subprocess
import sys
from collections import namedtuple
import meta_yank

ORDER_LIST = meta_yank.ORDER_LIST

META_YANK = '/home/jkb4y/work/results/meta/hg18/RegionYank/meta_yank.tbl'
EURMETA_YANK = '/home/jkb4y/work/results/eurmeta/hg18/RegionYank/eurmeta_yank.tbl'
OUTFOLDER = '/home/jkb4y/work/results/YankCompare/'
BFILE = '/home/jkb4y/work/data/2012Feb1/hg18/imm'

RegInfo = namedtuple('RegInfo',ORDER_LIST)
def read_yank(yank_loc):
    line1 = True
    yank_dict = dict()
    with open(yank_loc, mode="r") as yank:
        for line in yank:
            if line1:
                line1=False
                next
            else:
                line_split = line.strip().split()
                reg_info = RegInfo._make(line_split)
                yank_dict[reg_info.band]=reg_info
    return yank_dict









def plink_ld(bfile, plink_out,snp_one, snp_two):
    plink_cl_args = ['plink','--bfile', bfile,'--out', plink_out,
                     '--noweb', '--ld',snp_one,snp_two]
    p = subprocess.Popen(plink_cl_args, bufsize = 0, executable=None,stdin=None,
                         stdout=None,stderr=None, preexec_fn=None,close_fds=False,
                         shell=False,cwd=None,env=None, universal_newlines=False,
                         startupinfo=None, creationflags=0)
    p.wait()

def strip_name(yank_loc):
    head, tail = os.path.split(yank_loc)
    table_name = tail.replace('_yank.tbl','')
    return table_name


def read_log(log_loc):
    with open(log_loc, mode="r") as log:
        for line in log:
            if line.strip().startswith("R-sq"):
                print line
                line_list = line.strip().split()
                return line_list[2],line_list[5]
    return 'ERROR','ERROR'

def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global outfolder, annot, yank_one, yank_two
    global bfile
    
    yank_one = META_YANK
    yank_two = EURMETA_YANK
    outfolder = OUTFOLDER
    bfile = BFILE
    
    
    try: 
        opts, args = getopt.getopt(argv, "hb:o:",
                                   ["help","bfile=","outfolder=","y1=","yank-one=",
                                    "y2=","yank-two="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("--y1","--yank-one"):
            yank_one = arg
        elif opt in ("--y2","--yank-two"):
            yank_two = arg
        elif opt in ("-o","--outfolder"):
            outfolder = arg
        elif opt in ("-b","--bfile"):
            bfile = arg



def main(argv):
    global yank_one, yank_two, bfile, outfolder
    cl_arguments(argv)
    
    yank1_dict = read_yank(yank_one)
    yank2_dict = read_yank(yank_two)
    out_table = os.path.join(outfolder,'yank_comp.tbl')
    RegComp = namedtuple("RegComp",
                         "band,chro,start,end,sym,y1_lz,y1_im,y1_pos,y1_maf,y1_z,y1_p,y2_lz,y2_im,y2_pos,y2_maf,y2_z,y2_p,r2,D")
    comp_list = list()
    total = len(yank1_dict.keys())
    counter = 1
    print("There are {0} regions in {1}".format(total,yank_one))
    print("There are {0} regions in {1}".format(len(yank2_dict.keys()),yank_two))
    for band in yank1_dict.keys():
        print('''
        ****************************
Now contemplating region {0} of {1}.
        ****************************'''.format(counter,total))
        plink_out = os.path.join(outfolder,band)
        y1_snp = yank1_dict[band].lz
        y2_snp = yank2_dict[band].lz
        if not y1_snp == y2_snp:
            #plink_ld(bfile,plink_out,y1_snp,y2_snp)
            r2,d = read_log(plink_out+'.log')
        else:
            r2 = 'NA'
            d = 'NA'
        reg_comp = RegComp(band=band,chro=yank1_dict[band].chr,
                           start=yank1_dict[band].start,
                           end=yank1_dict[band].end,
                           sym=yank1_dict[band].sym,
                           y1_lz=yank1_dict[band].lz,
                           y1_im=yank1_dict[band].snp,
                           y1_pos=yank1_dict[band].pos,
                           y1_maf=yank1_dict[band].maf,
                           y1_z=yank1_dict[band].z,
                           y1_p=yank1_dict[band].p,
                           y2_lz=yank2_dict[band].lz,
                           y2_im=yank2_dict[band].snp,
                           y2_pos=yank2_dict[band].pos,
                           y2_maf=yank2_dict[band].maf,
                           y2_z=yank2_dict[band].z,
                           y2_p=yank2_dict[band].p,
                           r2=r2,
                           D=d)
        print reg_comp
        comp_list.append(reg_comp)
        counter = counter + 1
        
    yank1_name = strip_name(yank_one)
    yank2_name = strip_name(yank_two)
    title_list = ['chr_band','chr','start(Mb)','end(Mb)','gene_symbol',
                  yank1_name+'_snp',yank1_name+'_im_name',yank1_name+'_snp_pos',
                  yank1_name+'_maf',yank1_name+'_z_score',yank1_name+'_p_value',
                  yank2_name+'_snp',yank2_name+'_im_name',yank2_name+'_snp_pos',
                  yank2_name+'_maf',yank2_name+'_z_score',yank2_name+'_p_value',
                  'R-Sq',"D'"]
    with open(out_table, mode="w") as table:
        table.write('\t'.join(title_list)+'\n')
        for reg in comp_list:
            table.write('\t'.join(reg)+'\n')
            
        





if __name__ == '__main__':
    main(sys.argv[1:])
