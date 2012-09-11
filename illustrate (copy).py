#! ./
'''
@author: Jessica Bonnie
Created on February 7, 2012

'''
import getopt
import os
import sys
from collections import namedtuple

import pc_toolbox

import lz_meta
global table_loc
global outfolder
global region_loc
#global assoc
global population

TABLE_LOC = None
OUTFOLDER = None
OUT_FLAG = ''
CHROMOSOME = None
REGION_LOC = None
#DATA_SET = None
MAX_BP = None
POPULATION = None
#ASSOC = False
BUILD = None
CELL_TYPE = ''

CmdInfo = namedtuple("CmdInfo","flag,out,chro,max_bp,table,region,pop,assoc,fix")

def usage():
    print('''Usage goes here....        
USAGE: illustrate.py [FLAG] OBJECT
      FLAG                  DESCRIPTION                                     CURRENT DEFAULT
    -o, --outfolder         path of folder where results are to be written  {0}
    -c, --chromosome        chromosome number for which to draw macroview   {1}
    -t, --table             path to data table                              {2}
    -r, --region-list       path to region list                             {3}
    --max-bp                maximum number of bp contained in a macroview plot {4}
    --assoc
    --eqtl
    --pop                   LD population                                   {5}
    --build                 genome build: hg18 or hg19                      {6}
    -f, --flag              flag to add to output file names                {7}
    -h, --help              display this usage string
    '''.format(OUTFOLDER,CHROMOSOME,TABLE_LOC, REGION_LOC, lz_meta.MAX_BP,
               lz_meta.POPULATION,lz_meta.BUILD,
               OUT_FLAG))

def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global table_loc,chromosome, build
    global outfolder, out_flag, region_loc
    global hitspec
    #global assoc

    chromosome = None
    table_loc = TABLE_LOC
    outfolder = None
    out_flag = None
    region_loc = None
    #assoc = False
    #eqtl = False

    build = BUILD

    lz_cmdlist = list()
    
    try: 
        opts, args = getopt.getopt(argv, "ho:f:c:t:r:",
                                   ["help","outfolder=","flag=",
                                    "chromosome=","table=",
                                    "region-list=","pop=",
                                    "max-bp=","assoc","build=","hitspec",
                                    "chrom", "eqtl="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        cmd = [opt,arg]
        if opt in ("-h","--help"):
            usage()
            sys.exit()
##        elif opt in ("--pop","--max-bp","--build","--hitspec","--chrom"):
##            lz_cmdlist.extend(cmd)
##            print lz_cmdlist
        elif opt in ("-o","--outfolder"):
                outfolder = arg
                lz_cmdlist.extend(cmd)
        elif opt in ("-f","--flag"):
            out_flag = arg
            lz_cmdlist.extend(cmd)
        elif opt in ("-r","--region-list"):
            region_loc = arg
            lz_cmdlist.extend(cmd)
        elif opt in ("-c","--chromosome"):
            chromosome = arg
        elif opt in ("-t","--table"):
            table_loc = arg
            lz_cmdlist.extend(cmd)
        elif opt in ("--assoc"):
##            assoc = True
            lz_cmdlist.extend(["--assoc"])
        elif opt in ("--eqtl"):
            #eqtl = True
            lz_cmdlist.extend([opt, arg])
        elif opt in ("--pop","--max-bp","--build","--hitspec","--chrom"):
            lz_cmdlist.extend(cmd)
            print lz_cmdlist
    return lz_cmdlist


def create_command(lz_cmdlist, outfolder, logname, chromosome=None):
    cmd_list = ["python lz_meta.py"]
    cmd_list.extend(lz_cmdlist)
    if chromosome is not None:
        c_cmd = ["--chromosome",str(chromosome)]
        cmd_list.extend(c_cmd)
    #locate log folder
    lf = pc_toolbox.log_folder(outfolder)
    log_file = os.path.join(lf,logname)
    #add log command
    t_cmd = ["2>&1|tee",log_file]
    cmd_list.extend(t_cmd)
    
    cmd = ' '.join(cmd_list)
    print ('My log file is {0} and my command is {1}'.format(log_file,cmd))
    return cmd

def main():
    global chromosome, fix, table_loc, region_loc, out_flag, outfolder
    lz_cmdlist = cl_arguments(sys.argv[1:])
        
    print lz_cmdlist
    print outfolder
    #lz_cmdlist.extend(["--table",table_loc])
    #formulate name of logfile
    logname = ""
    of = ""
    if out_flag is not None:
        of = '_'+out_flag
                
    if chromosome is None and region_loc is None:
        for c in range(1,26):
            logname = 'chr'+str(c)+of+'.log'
            cmd = create_command(lz_cmdlist, outfolder, logname, c)
            os.system(cmd)
            sys.stdout.flush()
    elif chromosome is not None and region_loc is None:
        logname = 'chr'+chromosome+of+'.log'
        cmd = create_command(lz_cmdlist, outfolder, logname, chromosome)
        os.system(cmd)
        sys.stdout.flush()
    elif region_loc is not None:
        head, tail = os.path.split(region_loc)
        base, ext = os.path.splitext(tail)
        logname = "region_view"+of+'_'+base+".log"
        cmd = create_command(lz_cmdlist, outfolder, logname, chromosome)
        os.system(cmd)
        sys.stdout.flush()
    
if __name__ == '__main__':
    main()
