#! /usr/bin/python2.7
#! ./
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


REGION_LOC ='/home/jkb4y/work/Projects/UK/data/Region_Lists/T1D_regions_chr1.txt'
BP_FORM = 'mb'
CONDITION_LIST_FOLDER = None
#'/home/jkb4y/work/Projects/UK/data/Meta_Condition_Lists'
POSITION_FORM = 'mb'

CHRCOL = 0
STARTCOL = 2
ENDCOL = 3
GENECOL = 4
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
                            
    The following refer to the region-list:
        --chrcol        zero-based col # of chromosome column               {2}
        --genecol       zero-based col # of gene-symbol column              {3}
        --startcol      zero-based col # of position-start column           {4}
        --endcol        zero-based col # of position-end column             {5}
        --bp-form       start and end position units                        {6}
                        Options: bp, kb, mb

OPTIONS DIRECTLY PASSED TO PLINK_CONDITIONAL (and defaulted there)
    -o, --outfolder         path of folder where results are to be written  {7}
    -f, --flag              flag to add to output file names                {8}
    -p, --pbound            highest acceptable p-value (cut-off point)      {9}
    -l, --loop              maximum number of loops through PLINK           {10}
    -t, --test              which association test PLINK should run         {11}
    -s, --script            pathname of script to feed into PLINK           {12}
        --pheno             pathname of plink pheno file                    {13}
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
    ).format(REGION_LOC, CONDITION_LIST_FOLDER,str(CHRCOL),str(GENECOL),str(STARTCOL),
             str(ENDCOL),BP_FORM, plink_conditional.OUTPUT_FOLDER,
             plink_conditional.OUT_FLAG, str(plink_conditional.P_BOUND),
             str(plink_conditional.MAXLOOPS),plink_conditional.PLINK_TEST,
             plink_conditional.USER_SCRIPT_LOC,plink_conditional.PHENO_LOC))


def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global region_loc, chrcol, genecol, startcol, endcol
    global outflag, bp_form, condition_list_folder, pc_cmdlist
    
    pc_cmdlist = []
    
    region_loc = REGION_LOC
    chrcol = CHRCOL
    genecol = GENECOL
    startcol = STARTCOL
    endcol = ENDCOL
    outflag = OUTFLAG
    bp_form = BP_FORM
    condition_list_folder = CONDITION_LIST_FOLDER
    
    try: 
        opts, args = getopt.getopt(argv, "hf:r:l:p:t:s:o:",
                                   ["help","flag=","region-list=","loop=",
                                    "pbound=","test=","script=",
                                    "outfolder=",
                                    "chrcol=","startcol=",
                                    "endcol=","genecol=",
                                    "bp-form=","cfolder=",
                                    "pheno=","ci="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-l","--loop","-p","--pbound",
                     "-t","--test","-s","--script",
                     "--pheno","-o","--outfolder","--ci"):
            cmd = [opt,arg]
            sys.stdout.flush()
            pc_cmdlist.extend(cmd)
        elif opt in ("-f","--flag"):
            cmd = [opt, arg]
            sys.stdout.flush()
            pc_cmdlist.extend(cmd)
            outflag = arg
        elif opt in ("-r","--region-list"):
            region_loc = arg
        elif opt in ("--chrcol"):
            chrcol = int(arg)
        elif opt in ("--startcol"):
            startcol = int(arg)
        elif opt in ("--endcol"):
            endcol = int(arg)
        elif opt in ("--genecol"):
            genecol = int(arg)
        elif opt in ("--bp-form"):
            bp_form = arg.lower
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


def write_command(bp_form, region, pc_cmdlist):
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
    pos_start_flag = '--from-'+bp_form
    start_cmd = [pos_start_flag, start]
    cmd_list.extend(start_cmd)
    #include range end information
    pos_end_flag = '--to-'+bp_form
    print(pos_end_flag)
    end_cmd = [pos_end_flag, end]
    cmd_list.extend(end_cmd)
    #include refgene infromation
    rg_cmd = ['--refgene',gene]
    cmd_list.extend(rg_cmd)
    #include condition-list folder information
    if condition_list_folder is not None:
        cl_filename = 'Chr'+chromosome+'_'+gene+'.txt'
        cl_loc = os.path.join(condition_list_folder,cl_filename)
        cl_cmd = ["--condition-list",cl_loc]
        cmd_list.extend(cl_cmd)
        
    cmd = ' '.join(cmd_list)
    print("COMMAND TO PLINK CONDITIONAL IS: "+cmd)
    return cmd
    

def main():
    global region_loc, bp_form, pc_cmdlist
    #global chrcol, startcol, endcol, genecol
    print sys.argv
    cl_arguments(sys.argv[1:])
    print(pc_cmdlist)
##    gene_region_dict = pc_toolbox.create_region_dict(region_loc, POSITION_FORM)
##    key_list = gene_region_dict.keys()
##    print gene_region_dict
##    index_tup = (chrcol,startcol,endcol,genecol)
##    region_list = create_region_list(region_loc, index_tup)
    region_list = pc_toolbox.create_region_list(region_loc)
    for region in region_list:
        cmd = write_command(bp_form,region, pc_cmdlist)
        os.system(cmd)



if __name__=='__main__':
    main()




