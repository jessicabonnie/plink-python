#! /usr/bin/python2.6
#! ./
'''
Created on Jun 27, 2012

@author: Jessica Bonnie
'''
import os
import sys
import getopt
from collections import namedtuple

import pc_toolbox
import testsigs_workhorse

YANK_LOC = ''
OUTFOLDER = ''
CC_FOLDER = ''
YankInfo = namedtuple('YankInfo','start,end,chro,lz,pos,ID,p')
FREQ_LOC = '/home/jkb4y/work/data/2012Feb1/intersect_fam_uk/hg19/intersect_controls.frq'
LDFOLDER = '/home/jkb4y/cphgdesk_share/Projects/IMCHIP/Intersect_SNP_list/2012Feb1/eurmeta_LD/'
HIT_DICT_LOC = '/home/jkb4y/work/data/2012Feb1/intersect_fam_uk/TestSigSNPs/multihit.dict'

def read_yank(yank_loc):
    line1 = True
    with open(yank_loc, mode="r") as yank:
        yank_lines = list()
        for line in yank:
            if line1:
                index_dict = pc_toolbox.read_yank_titles(line)
                line1 = False
            else:
                lsplit = line.strip().split()
                if not lsplit == []:
                    yank_info = YankInfo(start=lsplit[index_dict['start']],
                                         end=lsplit[index_dict['end']],
                                         chro=lsplit[index_dict['chr']],
                                         lz=lsplit[index_dict['lz']],
                                         ID=lsplit[index_dict['id']],
                                         pos=lsplit[index_dict['pos']],
                                         p=lsplit[index_dict['p']])
                    yank_lines.append(yank_info)
    return yank_lines

def read_hit_dict(hit_dict_loc):
    hit_dict = dict()
    with open(hit_dict_loc) as hitty:
        for line in hitty:
            lsplit = line.strip().split()
            if not lsplit == []:
                hit_dict[lsplit[0]] = lsplit[1:]
    print hit_dict
    return hit_dict
                                     
def usage():
    print("usage goes here")   


def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global outfolder, ccfolder, yank_loc, single, ldfolder, freq_loc
    global out_flag, user_script_loc, build, horse_args, hit_dict_loc
    
    out_flag = ''
    yank_loc = None
    build = 'hg19'
    ccfolder = None
    outfolder = None
    user_script_loc = None
    horse_args = list()
    single = False
    ldfolder = LDFOLDER
    freq_loc = FREQ_LOC
    hit_dict_loc = HIT_DICT_LOC
    
    
    try: 
        opts, args = getopt.getopt(argv, "ho:s:f:",
                                   ["help","outfolder=","script=","flag=",
                                    "ccfolder=","yank=","build=","single",
                                    "freq=","ldfolder=","multihits="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-f","--flag"):
            out_flag = '_'+arg
        elif opt in ("-s","--script", "--build"):
            horse_args.extend((opt,arg))
            #user_script_loc = arg
        elif opt in ("-o","--outfolder"):
            outfolder = arg
        elif opt in ("multihits"):
            hit_dict_loc = arg
        elif opt in ("--ccfolder"):
            ccfolder = arg
        elif opt in ("--ldfolder"):
            ldfolder = arg
        elif opt in ("--freq"):
            freq_loc = arg
        elif opt in ("--single"):
            horse_args.append('--single')
            single = True
        
        elif opt in ("--yank"):
            yank_loc = arg
##        elif opt in ("--build"):
##            build = arg
##            print("BUILD IS SET TO: {0}!".format(arg))




def main(argv):
    global outfolder, ccfolder, yank_loc, ldfolder, freq_loc
    global out_flag, user_script_loc, build, horse_args, hit_dict_loc

    cl_arguments(argv)
    hit_dict = read_hit_dict(hit_dict_loc)
    yank_lines_list = read_yank(yank_loc)
##    region_out = '/home/jkb4y/work/data/Region_Lists/hg19/T1D_regions_hg19_05242012_sigsOnly.txt'
##    new_regionlist = open(region_out, mode = "w")
##    new_regionlist.write('\t'.join(['gene_chr','Chr_band',
##                          'region_start','region_end',
##                          'gene_symbol','band_title','region_id'])+'\n')
    for info in yank_lines_list:
        if not float(info.p) > 3.5e-7:
            chromosome = info.chro
            start = info.start
            end = info.end
            lz = info.lz
            ID = info.ID
            band_title = info.ID.replace('_z','').replace('_a','').replace('_b','')
            
##            new_regionlist.write('\t'.join([chromosome, ID, start,
##                                            end, '.',band_title,ID])+'\n')
            
            assoc_name = ID + '_~0_SNPsOut.assoc.logistic'
            assoc_loc = os.path.join(ccfolder, 'chr{0}'.format(chromosome),assoc_name)
        
            cmd_list = ['--region-id',ID,'--chromosome',chromosome,
                        '--start',start,'--end',end,
                        #'--snpstar',lz,
                        '--assoc',assoc_loc]#'--build', build,
                        #'--script',user_script_loc,'--outfolder',outfolder]
            if single:
                if not 'Singles' in outfolder:
                    print("in if statement")
                    new_out = os.path.join(outfolder, 'Singles')
                    if not os.path.exists(new_out):
                        os.makedirs(new_out)
                    outfolder = new_out
            cmd_list.extend(['--outfolder',outfolder])
            cmd_list.extend(['--ldfolder',ldfolder,'--freq',freq_loc])
            cmd_list.extend(horse_args)
            hit_index = 1
            if ID in hit_dict:
                for snp in hit_dict[ID]:
                    cmd_list_prime = []
                    cmd_list_prime.extend(cmd_list)
                    if hit_index == 1:
                        cmd_list_prime.extend(['--hit1',','.join(hit_dict[ID])])
                    cmd_list_prime.extend(['--snpstar',snp,'--multi', hit_index])
                    print cmd_list_prime
                    testsigs_workhorse.main(cmd_list_prime)
                    hit_index = hit_index + 1
                    
            else:
                cmd_list.extend(['--snpstar',lz])
                print cmd_list
                testsigs_workhorse.main(cmd_list)
                
##    new_regionlist.close()
    


if __name__ == '__main__':
    main(sys.argv[1:])
