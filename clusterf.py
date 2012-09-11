#! /usr/bin/python2.6
#! ./
'''
@author: Jessica Bonnie
Created on Mar 1, 2011

'''
import os
import sys
import getopt

import plink_conditional
import pc_toolbox

global region_loc
global condition_list_folder

BUILD = 'hg18'
REGION_LOC = None
#'/home/jkb4y/work/data/Region_Lists/hg18/Chromosomes.txt'
CONDITION_LIST_FOLDER = None
OUTFLAG = None

def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global region_loc, cluster, cflag
    global outflag, condition_list_folder, pc_cmdlist, build
    
    pc_cmdlist = []
    
    region_loc = REGION_LOC
    outflag = OUTFLAG
    condition_list_folder = CONDITION_LIST_FOLDER
    build = BUILD
    cluster = False
    cflag = ''
    try: 
        opts, args = getopt.getopt(argv, "hf:r:l:p:t:s:o:",
                                   ["help","flag=","region-list=","loop=",
                                    "pbound=","test=","script=",
                                    "outfolder=",
                                    "cfolder=",
                                    "pheno=","ci=","build=",
                                    "cluster"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-l","--loop","-p","--pbound",
                     "-t","--test","-s","--script",
                     "--pheno","-o","--outfolder","--ci","--build"):
            cmd = [opt,arg]
            sys.stdout.flush()
            pc_cmdlist.extend(cmd)
        elif opt in ("-f","--flag"):
            cmd = [opt, arg]
            sys.stdout.flush()
            pc_cmdlist.extend(cmd)
            outflag = arg
            cflag + '_'+arg
        elif opt in ("-r","--region-list"):
            region_loc = arg
        elif opt in ("--cluster"):
            cluster = True
##        elif opt in ("--endcol"):
##            endcol = int(arg)
##        elif opt in ("--genecol"):
##            genecol = int(arg)
##        elif opt in ("--bp-form"):
##            bp_form = arg.lower
        elif opt in ("--cfolder"):
            condition_list_folder = arg
            head, tail = os.path.split(condition_list_folder)
            cflag = cflag + '_'+tail.replace('condition_lists','').replace('lz','').replace('_','')
            print cflag
            
            
            

##def create_region_col_tup():
##    global chrcol, startcol, endcol, genecol
##    col_index_tup = (chrcol,startcol,endcol,genecol)
##    return col_index_tup


##def create_region_list(region_loc, region_col_indices):
##    region_list = []
##    with open(region_loc, mode = "r") as region_file:
##        region_file.next()
##        for line in region_file:
##            line_split = line.strip().split()
##            chromosome = str(line_split[region_col_indices[0]])
##            start_str = str(line_split[region_col_indices[1]])
##            end_str = str(line_split[region_col_indices[2]])
##            RegInfo = namedtuple('RegInfo', 'chro,start,end,sym')
##            gene_sym = pc_toolbox.correct_gene_sym(line_split[region_col_indices[3]],chromosome)
##            region_tup = RegInfo(chro=chromosome,start=start_str,end=end_str,sym=gene_sym)
##            region_list.append(region_tup)
##    return region_list


def write_command(region, pc_cmdlist):
    global outflag, condition_list_folder
    
    chromosome = region.chro
    start = region.start
    end = region.end
    gene = region.sym
    band = region.band
    #start with a call to plink_conditional.py
    cmd_list = ["python plink_conditional.py"]
    cmd_list.extend(pc_cmdlist)
    #include chromosome information
    chr_cmd = ["--chromosome",chromosome]
    cmd_list.extend(chr_cmd)
    #include range start information
    pos_start_flag = '--from-mb'
    start_cmd = [pos_start_flag, start]
    cmd_list.extend(start_cmd)
    #include range end information
    pos_end_flag = '--to-mb'
    end_cmd = [pos_end_flag, end]
    cmd_list.extend(end_cmd)
    #include refgene information
    rg_cmd = ['--refgene',gene]
    cmd_list.extend(rg_cmd)
    #include chrband information
    cb_cmd = ['--chrband',band]
    cmd_list.extend(cb_cmd)
    #include condition-list folder information
    if condition_list_folder is not None:
         cl_filename = band + '.txt'
##       cl_filename = 'Chr'+chromosome+'_'+gene+'.txt'
         cl_loc = os.path.join(condition_list_folder,cl_filename)
         cl_cmd = ["--condition-list",cl_loc]
         cmd_list.extend(cl_cmd)
        
    cmd = ' '.join(cmd_list)
    print("COMMAND TO PLINK CONDITIONAL IS: "+cmd)
    return cmd
    

def main(argv):
     global region_loc, pc_cmdlist, cluster, cflag
     print sys.argv
     cl_arguments(argv)
     print(pc_cmdlist)
##    gene_region_dict = pc_toolbox.create_region_dict(region_loc, POSITION_FORM)
##    key_list = gene_region_dict.keys()
##    print gene_region_dict
##    index_tup = (chrcol,startcol,endcol,genecol)
##    region_list = create_region_list(region_loc, index_tup)
     region_list = pc_toolbox.create_region_list(region_loc)
     for region in region_list:
          cmd = write_command(region, pc_cmdlist)
          if cluster:
               bash_loc = './PBS_scripts/'+region+'_'+cflag+'.sh'
               print bash_loc
               #write_bash(bash_loc, cmd)
               #os.system('qsub {0}.'.format(bash_loc))
          else:
               os.system(cmd)



if __name__=='__main__':
    main(sys.argv[1:])





