#! /usr/bin/python2.7
#! /usr/bin/python2.6
#! ./
#PBS -l select=1:ncpus=3:mem=1gb
'''
@author: Jessica Bonnie
Created on Nov 3, 2011

'''
import os
import sys
import getopt
import plink_conditional
import pc_toolbox
from collections import namedtuple


global region_loc
global bp_form
global condition_list_folder
global interrupt

BUILD = 'hg19'
REGION_LOC = None
#'/home/jkb4y/work/data/Region_Lists/hg18/Chromosomes.txt'
CONDITION_LIST_FOLDER = None
OUTFLAG = None


##
##def correct_gene_symbol(gene_symbol, chromosome):
##    ''' Determine if reference gene symbol is appropriate for dictionary use,
##        and change where appropriate.
##        Args:
##            gene_symbol -- gene_symbol listed in the region list
##            chromosome -- chromosome listed in region list
##        Returns:
##            new_symbol -- symbol to be used in dictionary
##    
##    '''
##    new_symbol = gene_symbol
##    if gene_symbol=='0' or gene_symbol =='No' or gene_symbol == 'no' or gene_symbol == 'NA':
##        new_symbol = 'nogene-'+ chromosome
##    elif gene_symbol=='multiple' or gene_symbol=='Multiple':
##        new_symbol = 'multiple-'+chromosome
##    return new_symbol


def usage():
     print(

    ('''
USAGE: region_wrap.py [FLAG] OBJECT
      FLAG                  DESCRIPTION                                     CURRENT DEFAULT
    -r, --region-list       path location of region-list                    {0}
    --cfolder               path location of folder of
                                condition-lists produced by meta_yank.py    {1}

OPTIONS DIRECTLY PASSED TO PLINK_CONDITIONAL (and defaulted there)
    -o, --outfolder         path of folder where results are to be written  {2}
    -f, --flag              flag to add to output file names                {3}
    -p, --pbound            highest acceptable p-value (cut-off point)      {4}
    -l, --loop              maximum number of loops through PLINK           {5}
    -t, --test              which association test PLINK should run         {5}
    -s, --script            pathname of script to feed into PLINK           {7}
        --pheno             pathname of plink pheno file                    {8}
    -h, --help              display this usage string
    
OBJECT Options:
    FLAG            OPTION          DESCRIPTION
  [-t] ----->       assoc           Case/Control for QTL association
                    fisher          Fisher's exact (allelic) test
                    model           Cochran-Armitage and full-model C/C association
                    mh              Cochran-Mantel-Haenszel SNPxDISEASE|STRATA
                    linear          Test for quantitative traits and multiple covariates
                    logistic        Test for disease traits and multiple covariates         
    '''
    ).format(REGION_LOC, CONDITION_LIST_FOLDER,
             plink_conditional.OUTPUT_FOLDER,
             plink_conditional.OUT_FLAG, str(plink_conditional.P_BOUND),
             str(plink_conditional.MAXLOOPS),plink_conditional.PLINK_TEST,
             plink_conditional.USER_SCRIPT_LOC,plink_conditional.PHENO_LOC))


def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    
    global region_loc, interrupt
    global outflag, condition_list_folder, pc_cmdlist, build

    
    pc_cmdlist = []
    
    region_loc = REGION_LOC
    outflag = OUTFLAG
    condition_list_folder = CONDITION_LIST_FOLDER
    build = BUILD
    interrupt = None
    
    try:
          opts, args = getopt.getopt(argv, "hf:r:l:p:t:s:o:",
                                     ["help","flag=","region-list=","loop=",
                                      "pbound=","test=","script=","outfolder=",
                                      "cfolder=","pheno=","ci=",
                                      "build=","interrupt="])
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
        elif opt in ("--interrupt"):
            sys.stdout.flush()
            pc_cmdlist.extend([opt])
            interrupt = arg
        elif opt in ("-r","--region-list"):
            region_loc = arg
        elif opt in ("--cfolder"):
            condition_list_folder = arg
            
            

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
    global outflag, condition_list_folder, interrupt
    
    chromosome = region.chro
    start = region.start
    end = region.end
    ID = region.ID
    title = region.title
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
##    rg_cmd = ['--refgene',gene]
##    cmd_list.extend(rg_cmd)
    #include chrband information
    cb_cmd = ['--chrband',ID]
    cmd_list.extend(cb_cmd)
    #include condition-list folder information
    if interrupt is not None:
         cl = interrupt
         cl_cmd = ["--condition-list",cl]
         cmd_list.extend(cl_cmd)
    elif condition_list_folder is not None:
         cl_filename = ID + '.txt'
##       cl_filename = 'Chr'+chromosome+'_'+gene+'.txt'
         cl_loc = os.path.join(condition_list_folder,cl_filename)
         cl_cmd = ["--condition-list",cl_loc]
         cmd_list.extend(cl_cmd)
        
    cmd = ' '.join(cmd_list)
    print("COMMAND TO PLINK CONDITIONAL IS: "+cmd)
    return cmd
    

def main(argv):
    global region_loc, pc_cmdlist
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
         os.system(cmd)



if __name__=='__main__':
    main(sys.argv[1:])




