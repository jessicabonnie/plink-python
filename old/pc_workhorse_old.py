#! /usr/bin/python2.7
#! ./
'''
Created on Sep 14, 2011

@author: Jessica Bonnie
'''

import os
import subprocess
import getopt
import sys
import time
import plink_conditional
import shutil
import re
import math
from collections import namedtuple

global loopcount	#how many times looped thru plink
global SNP_loc		#file location of list of SNPs to ignore

global user_script_loc  #location of user's supplied script

global output_folder    #user's supplied output folder
global out_flag         #user's supplied flag for naming files


global snp_col		#index location of SNPs column in assoc file from PLINK
global p_col		#index location of P column in assoc file from PLINK
global firstloop        #boolean indicating that this is the first time thru, even if loopcount is not 0

global p_bound          #user's supplied lowest P-value of interest
global loopagain        #boolean indicating whether another loop will follow current one
global maxloops         #user's supplied maximum number of times to loop through plink
global plink_test       #user's supplied choice of plink association test
global range_start_bp   #translation of user's supplied start to range on chromosome
global range_end_bp     #translation of user's supplied end to range on chromosome
global refgene          #user's supplied reference gene name (used in naming)
global condition_list   #user's supplied location of condition-list
global chromosome       #user's supplied chromosome number

global pheno_tag        #current phenotype on which pc_workhorse is being run




#DEFAULT VALUES
#N.B. > THESE ARE NOT TO BE TOUCHED BY USER. IF USER WISHES TO ADAPT
#DEFAULTS, THEY MAY BE ADAPTED WITHIN THE WRAPPER: plink_conditional.py

OUT_FLAG = ''
REFGENE = None
CONDITION_LIST = None
PHENO_TAG = ''



def identify_cols(first_line_list, plink_test, c_interval):
    '''assigns locations of P and SNP columns to global variables

    Keyword arguments:
    first_line_list -- list created from the first line of data file,
        containing column titles
    Returns - tuple containing the indices of the P and SNP column titles

    '''
    
    p_index = first_line_list.index('P')
    snp_index = first_line_list.index('SNP')
    a1_index = first_line_list.index('A1')
    or_index = None
    stat_index = None
    ci_hi_index = None
    ci_lo_index = None
    if plink_test in ('logistic', 'fisher','assoc', 'linear'):
        or_index = first_line_list.index('OR')
    if plink_test in ('linear', 'assoc'):
        or_index = first_line_list.index('BETA')
    if plink_test in ('logistic','linear'):
        stat_index = first_line_list.index('STAT')
    if c_interval is not None:
        ci = str(int(float(c_interval)*100))
        ci_lo_index = first_line_list.index('L'+ci)
        ci_hi_index = first_line_list.index('U'+ci)
    index_dict = {'p':p_index,
                  'snp':snp_index,
                  'or':or_index,
                  't':stat_index,
                  'hi':ci_lo_index,
                  'lo':ci_hi_index,
                  'a1':a1_index}
    return index_dict

def find_min(tuple_list):
    from operator import itemgetter
    
    sorted_tlist = sorted(tuple_list, key=itemgetter(0,2))
    


def findSNPWithSmallestPValue(tuple_list):
    '''returns the SNP with the lowest p-value

    Keyword arguments:
    tuple_list -- list of tuples of the form (p-value,SNP)
    Returns - SNP with the lowest p-value

    '''
    global loopcount, condition_list, p_bound, loopagain
    tuplet = min(tuple_list)
    find_min(tuple_list)
    if tuplet.p > p_bound:
        loopagain = False
    cl_indicator = ''
    if condition_list is not None:
        cl_indicator = '*'
    stat_label = ''
    or_label = ''
    if not tuplet[3] == '':
        if plink_test in ('logistic', 'linear'):
            stat_label = 'Coefficient T-Statistic: {0}'.format(str(tuplet[3]))
    if not tuplet[2] == '':
        if plink_test in ('linear'):
            or_label = 'Regression Coefficient: {0}'.format(str(tuplet[2]))
        if plink_test in ('logistic','fisher','assoc'):
            or_label = 'Odds Ratio: {0}'.format(str(tuplet[2]))
        or_label = or_label + '\nConfidence Interval: {0} - {1}'.format(str(tuplet[4]),
                                                                         str(tuplet[5]))
        
    print('''
*****************************************************
                In loop {0}
    {1} has the lowest P value --------> {2}
%%%         {1}     {2}     {0}
{3}
{4}
{5}
*****************************************************

'''.format(str(loopcount)+cl_indicator,str(tuplet[1]),str(tuplet[0]),
           or_label, stat_label, 'A1: '+ str(tuplet.a1)))
    snp_addition = tuplet.snp
    if not loopagain:
	print('''
**************************************************************
    **************************************************
    There are no more SNPs with P-Values less than {0}!
        This was the last loop through PLINK!
:::                     {1}
    **************************************************
***************************************************************
'''.format(str(p_bound),str(loopcount)+cl_indicator))
	snp_addition = ''
    sys.stdout.flush()
    return snp_addition

def tuple_helper(line_list, index):
    #print line_list
    #print index
    item = ''
    if index is not None:
        item = line_list[int(index)]
    #print item
    sys.stdout.flush()
    return item
    
def get_tuple(line_list, index_dict):
    '''creates named tuple of the form ('p'=p-value,
                                        'snp'=SNP,
                                        'or'=odds ratio,
                                        't'=t-statistic,
                                        'lo'=lower confidence interval edge,
                                        'hi' = upper ci edge,
                                        'a1' = A1 allele)

    Keyword arguments:
    line_list -- line split into its list of elements

    Returns - tuple of the form (p-value, SNP, odds ratio, t-statistic)

    '''
    p_store = line_list[index_dict['p']]
    if is_number(line_list[index_dict['p']] ):
        p_store = float(line_list[index_dict['p']])
    LineInfo = namedtuple('LineInfo', 'p,snp,OR,t,lo,hi,a1')
    log_tuple = LineInfo(p=p_store,
                         snp=line_list[index_dict['snp']],
                         OR=tuple_helper(line_list,index_dict['or']),
                         t=tuple_helper(line_list,index_dict['t']),
                         lo=tuple_helper(line_list,index_dict['lo']),
                         hi=tuple_helper(line_list,index_dict['hi']),
                         a1=tuple_helper(line_list,index_dict['a1']))
    return log_tuple
    

def plink(script_loc):
    '''runs a plink's logistic association test

    Keyword arguments:
    script_loc -- location of plink script

    Returns - plink subprocess
    '''
    p=subprocess.Popen(('plink','--script',script_loc),
                           bufsize = 0, executable=None,stdin=None,stdout=None,
                           stderr=None, preexec_fn=None,close_fds=False,shell=False,
                           cwd=None,env=None, universal_newlines=False,startupinfo=None,
                           creationflags=0)
    p.wait()

    

def output_subfolder(output_folder, chromosome):
    '''Creates an output subfolder based on the output folder and the chromosome number.
        Args:
            output_folder -- filepath of output mother folder
            chromosome -- chromosome number in string form
    '''
    output_subfolder = os.path.join(output_folder, 'chr'+chromosome)
    if not os.path.exists(output_subfolder):
        os.makedirs(output_subfolder)
    return output_subfolder


def is_number(s):
    '''checks a string to determine if it is a number

    '''
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def read_assoc(assoc_loc, SNP_loc,c_interval,plink_test, exp):
    global p_bound, index_dict, firstloop
    tuple_list = []
    with open(assoc_loc, mode = 'r') as assoc_file:
        line1_check = True                                        
        for assoc_line in assoc_file:                           #open the assoc file to be read
            assoc_list = assoc_line.split()
            
            if line1_check:
                if firstloop:
                    index_dict = identify_cols(assoc_list, plink_test, c_interval)
                    firstloop = False
                if loopcount == 0:
                    if os.path.exists(SNP_loc):
                        clear = open(SNP_loc, mode='w')
                        clear.close()
                line1_check = False
            else:
                tuplet = get_tuple(assoc_list, index_dict)
                if is_number(tuplet.p) and tuplet.p < float(p_bound)* math.pow(10,exp):
                        tuple_list.append(tuplet)
    return tuple_list

def find_and_record(assoc_loc, SNP_loc, condition_list,plink_test,c_interval):
    '''locates the SNP with the lowest p-value and appends it to a file

    Keyword arguments:
    assoc_loc -- data file from which the p-values and SNPs are to be read
    GLOBALS utilized: SNP_loc, p_bound, loopagain, loopcount, firstloop

    Returns - file with list of lowest p-valued snps

    '''
    global p_bound,  loopcount, firstloop, index_dict
    
    exp = 2
    cl_indicator = ''
    if condition_list is not None:
        cl_indicator = '*'
    tuple_list = read_assoc(assoc_loc, SNP_loc, c_interval, plink_test, exp)
##    with open(assoc_loc, mode = 'r') as assoc_file:
##        line1_check = True                                        
##        for assoc_line in assoc_file:                           #open the assoc file to be read
##            assoc_list = assoc_line.split()
##            
##            if line1_check:
##                if firstloop:
##                    index_dict = identify_cols(assoc_list, plink_test, c_interval)
##                    firstloop = False
##                if loopcount == 0:
##                    if os.path.exists(SNP_loc):
##                        clear = open(SNP_loc, mode='w')
##                        clear.close()
##                line1_check = False
##            else:
##                tuplet = get_tuple(assoc_list, index_dict)
##                if is_number(tuplet.p) and tuplet.p < float(p_bound)*100:
##                        tuple_list.append(tuplet)
                        
    while len(tuple_list) == 0:
        exp = exp + 2
        tuple_list = read_assoc(assoc_loc, SNP_loc, c_interval, plink_test, exp)
        
##        loopagain = False
##        print('''
##**************************************************************
##    **************************************************
##    There are no more SNPs with P-Values less than {0}!
##        This was the last loop through PLINK!
##:::                     {1}
##    **************************************************
##***************************************************************
##'''.format(str(p_bound),str(loopcount)+cl_indicator))
##        sys.stdout.flush()
##        return
    with open(SNP_loc, mode = 'a') as SNP_file:     #open/create the SNP file
            SNP_file.write(findSNPWithSmallestPValue(tuple_list))
            SNP_file.write('\n')

def usage():
    print('''
ERROR.  DO NOT ACCESS PC_WORKHORSE.PY DIRECTLY.
        USE WRAPPER FUNCTION: PC_CONDITIONAL.PY or REGION_WRAP.PY
''')
##USAGE: pc_workhorse.py [FLAG] OBJECT
##      FLAG                  DESCRIPTION                                     CURRENT DEFAULT
##    -o, --outfolder         path of folder where results are to be written  {0}
##    -f, --flag              flag to add to output file names                {1}
##    -c, --chromosome        number of chromosome on which to run test       {2}
##        --from-bp	    start of range on chromosome in bp		    
##	--to-bp		    end of range on chromosome in bp		    
##	--from-kb           start of range on chromosome in kb
##        --to-kb             end of range on chromosome in kb
##        --from-mb           start of range on chromosome in mb
##        --to-mb             end of range on chromosome in mb
##	--refgene	    Gene Name of Interest (for naming purposes)	    
##    -p, --pbound            highest acceptable p-value (cut-off point)      {3}
##    -l, --loop              maximum number of loops through PLINK           {4}
##    -t, --test              which association test PLINK should run         {5}
##    -s, --script            pathname of script to feed into PLINK           {6}
##        --pheno-tag         Name of pheno-type of interest in current run
##        --condition-list    pathname of previously compiled condition-list
##    -h, --help              display this usage string
##    
##OBJECT Options:
##    FLAG            OPTION          DESCRIPTION
##  [-t] ----->       assoc           Case/Control for QTL association
##                    fisher          Fisher's exact (allelic) test
##                    model           Cochran-Armitage and full-model C/C association
##                    mh              Cochran-Mantel-Haenszel SNPxDISEASE|STRATA
##                    linear          Test for quantitative traits and multiple covariates
##                    logistic        Test for disease traits and multiple covariates         
##    '''.format(plink_conditional.OUTPUT_FOLDER,OUT_FLAG,plink_conditional.CHROMOSOME,
##             str(plink_conditional.P_BOUND),plink_conditional.MAXLOOPS,
##             plink_conditional.PLINK_TEST,plink_conditional.USER_SCRIPT_LOC))

def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global output_folder, plink_test, p_bound, chromosome
    global maxloops, out_flag, user_script_loc
    global range_start_bp, range_end_bp, refgene
    global condition_list, pheno_tag, c_interval
    
    out_flag = OUT_FLAG
    refgene = REFGENE
    condition_list = CONDITION_LIST
    pheno_tag = PHENO_TAG
    c_interval = None
    plink_test = None
    
    try: 
        opts, args = getopt.getopt(argv, "ho:s:f:t:c:p:l:g:",
                                   ["help","outfolder=","script=","flag=","test=",
                                    "chromosome=","pbound=","loop=","refgene=",
                                    "from-bp=","to-bp=",
                                    "from-kb=","to-kb=","from-mb=","to-mb=",
                                    "condition-list=","pheno-tag=","ci="])
    except getopt.GetoptError:
        plink_conditional.usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-f","--flag"):
            out_flag = '_'+arg
        elif opt in ("-s","--script"):
            user_script_loc = arg
        elif opt in ("-g","--refgene"):
            refgene = arg
        elif opt in ("-o","--outfolder"):
            output_folder = arg
        elif opt in ("-t","--test"):
            plink_test = arg
        elif opt in ("-c","--chromosome"):
            chromosome = arg
        elif opt in ("--from-bp"):
            range_start_bp = int(arg)
        elif opt in ("--to-bp"):
            range_end_bp = int(arg)
        elif opt in ("--from-kb"):
            range_start_bp = int(1000 * float(arg))
        elif opt in ("--to-kb"):
            range_end_bp = int(1000 * float(arg))
        elif opt in ("--from-mb"):
            range_start_bp = int(1e6 * float(arg))
        elif opt in ("--to-mb"):
            range_end_bp = int(1e6 * float(arg))
        elif opt in ("-p","--pbound"):
            p_bound = float(arg) 
        elif opt in ("-l","--loop"):
            maxloops = int(arg)
        elif opt in ("--condition-list"):
            condition_list = arg
        elif opt in ("--pheno-tag"):
            pheno_tag = arg
        elif opt in ("--ci"):
            c_interval = arg


def get_user_script_text(user_script_loc):
    '''Obtain contents of user script.
        Args:
            user_script_loc -- filepath of user script
        Returns:
            uss -- user script string
    '''
    uss = ''
    with open(user_script_loc, mode='r') as user_script:
        for line in user_script:
            uss = uss + line + '\n'
    return uss

def write_script(output_file, new_script_loc):
    global SNP_loc, loopcount,chromosome, range_start_bp
    global range_end_bp, user_script_loc, c_interval
    user_script_string = get_user_script_text(user_script_loc)
    with open(new_script_loc, mode="w") as new_script:
        new_script.write(user_script_string + '\n')
        new_script.write('--' + plink_test +'\n' +
                         '--chr ' + chromosome +'\n'+
                         '--from-bp ' + str(range_start_bp) +'\n'+
                         '--to-bp ' + str(range_end_bp) +'\n'+
                         '--out ' + output_file)
        if not loopcount == 0:
            new_script.write('\n--condition-list '+SNP_loc)
        if c_interval is not None:
            new_script.write('\n--ci '+ c_interval)
##    if loopcount == 0:
##        with open(new_script_loc,mode='w') as new_script:
##            new_script.write(user_script_string +'\n' +
##                             '--' + plink_test +'\n' +
##                             '--chr ' + chromosome +'\n'+
##                             '--from-bp ' + str(range_start_bp) +'\n'+
##                             '--to-bp ' + str(range_end_bp) +'\n'+
##                             '--out ' + output_file)
##            
##    else:
##        with open(new_script_loc,mode='w') as new_script:
##            new_script.write(user_script_string +'\n' +
##                             '--' + plink_test +'\n' +
##                             '--chr ' + chromosome +'\n'+
##                             '--from-bp ' + str(range_start_bp) +'\n'+
##                             '--to-bp ' + str(range_end_bp) +'\n'+
##                             '--out ' + output_file + '\n'
##                             '--condition-list ' + SNP_loc )

def convert_bp_to_mb(bp):
    mb = str(round((float(bp)*1e-6),2))
    return mb
            

def range_label(chromosome, refgene, range_start_bp,range_end_bp):
    mb_start = convert_bp_to_mb(range_start_bp)
    mb_end = convert_bp_to_mb(range_end_bp)
    if refgene is None:
        rlabel = 'Chr' + chromosome + ':'+mb_start+'-'+mb_end
    else:
        rlabel = 'Chr' + chromosome + ':' + refgene
    return rlabel

def locuszoom(assoc_loc, range_label, chromosome,range_start_bp,range_end_bp,pheno_tag, label):
    global loopcount
##    global refgene, out_flag, condition_list
    assoc_head, assoc_tail = os.path.split(assoc_loc)
    title = '{0} Region Case/Control: {1} SNPs Conditioned Out'.format(range_label,str(loopcount))
    if not pheno_tag == '':
        title = title + '\n Conditioned on Expression of {0}'.format(pheno_tag)
    if '~' in label:
        title = title + '\n(Condition List Provided at Start)'
        
    lz_cl_args = ['locuszoom','--metal', str(assoc_loc),
                  '--chr',str(chromosome),'--start', str(range_start_bp),
                  '--end', str(range_end_bp), '--pvalcol','P',
                  '--markercol','SNP','--delim','whitespace', '--prefix',
		  label,'--plotonly','--verbose','--no-date',
                  '--snpset','None','showRecomb=FALSE','title='+title]
    p = subprocess.Popen(lz_cl_args, bufsize = 0, executable=None,stdin=None,
                         stdout=None,stderr=None, preexec_fn=None,close_fds=False,
                         shell=False,cwd=None,env=None, universal_newlines=False,
                         startupinfo=None, creationflags=0)
    p.wait()



def determine_extension(assoc_test):
    extension = ''
    if assoc_test == 'logistic':
        extension = extension+'.assoc.logistic'
    elif assoc_test == 'assoc':
        extension = extension+'.assoc'
    elif assoc_test == 'linear':
        extension = extension+'.assoc.linear'
    elif assoc_test == 'hap-assoc':
        extension = extension+'.assoc.hap'
    elif assoc_test == 'model':
        extension = extension+'.model'
    elif assoc_test == 'fisher':
        extension = extension+'.fisher'
    elif assoc_test == 'mh':
        extension = extension+'.cmh'
    else:
        print("Test produces unknown file extension, please adapt code.")
        exit(2)
    return extension

def print_snp_list(SNP_loc):
    global loopcount
    if not loopcount == 0:
        print('''
**********************************************************************
    ************************************************************
    The data will be conditioned on the following SNPs:
''')
        with open(SNP_loc, mode = "r") as snippy:
            for snp in snippy:
                print('&&&\t\t'+snp)
        print('''
        ******************************************************
**********************************************************************
''')
    sys.stdout.flush()

def given_condition_list(SNP_loc, condition_list):
    global loopcount
    print('''
**********************************************************************
    ************************************************************
    NOTE: CONDITION LIST has been provided.
    ************************************************************
**********************************************************************
''')
    sys.stdout.flush()
    counter = 0
    list_loc = condition_list
    if condition_list == SNP_loc:
        (basepath, ext) = os.path.splitext(SNP_loc)
        new_loc = os.path.join(basepath + '~input'+ext)
        print('''
The provided condition-list path, {0}, is the same as the output condition-list,
the contents of the provided list will be moved to {1}, and will also be included
in the output list.
'''.format(condition_list, new_loc))
        sys.stdout.flush()
        shutil.copy(SNP_loc, new_loc)
        list_loc = new_loc
    print('''
**********************************************************************
        ************************************************************
        The following SNPs were found in the condition-list:
''')
    with open(SNP_loc, mode = 'w') as SNP_file:
        with open(list_loc, mode = 'r') as condition_file:
            for condition in condition_file:
                if not condition.split() == []:
                    SNP_file.write(condition + '\n')
                    print('\t\t'+condition)
                    counter = counter + 1
    loopcount = counter
    print('''
        Number of SNPs in original condition-list: {0}
        ******************************************************
**********************************************************************
'''.format(loopcount))
    sys.stdout.flush()

def main():
    global loopcount, output_folder, plink_test, chromosome, refgene
    global maxloops, loopagain, SNP_loc, range_start_bp,range_end_bp
    global out_flag, condition_list, firstloop, pheno_tag, c_interval

    firstloop = True
    loopagain = True
    loopcount = 0
    cl_arguments(sys.argv[1:])
    
    label = ''
    mb_start = convert_bp_to_mb(range_start_bp)
    mb_end = convert_bp_to_mb(range_end_bp)
    cl_indicator = ''
    rgene = refgene
    ptag = ''
    if not pheno_tag == '':
        ptag = '_'+pheno_tag

    rangelabel = range_label(chromosome, refgene, range_start_bp,range_end_bp)
    if condition_list is not None:
        cl_indicator = '~'
    if refgene is None:
        rgene = 'NA'
    label = rangelabel + ptag+out_flag +'_'+cl_indicator
    of = '--'
    if not out_flag == '':
        of_pretty_list = out_flag.split('_')
        of = of_pretty_list[1]
    sys.stdout.write('''

******************************************************************************
        ***********************************************************
                    Beginning pc_workhorse.py on {0} at
                                    {1}
$$$     {2}             {3}             {4}         {5}         {6}
        ***********************************************************
******************************************************************************

'''.format(label, time.strftime("%a,%c"),chromosome,rgene,mb_start,mb_end, of+ptag))
    sys.stdout.flush()
    output_folder = str(output_subfolder(output_folder,chromosome))
    os.chdir(output_folder)
    assoc_extension = determine_extension(plink_test)
    
    SNP_loc=os.path.join(output_folder, label + 'leastP_SNPs.txt')
    if condition_list is not None:
        given_condition_list(SNP_loc, condition_list)
        
    
    while (loopcount < maxloops and loopagain==True):
        print_snp_list(SNP_loc)
        assoc_file_basename = label + str(loopcount)+'_SNPsOut'
        assoc_file_name = assoc_file_basename + assoc_extension
        script_name = label+ 'script_' + str(loopcount) + '.txt'
        write_script(assoc_file_basename, script_name)
        plink(script_name)
        find_and_record(assoc_file_name, SNP_loc, condition_list, plink_test, c_interval)
        sys.stdout.flush()
        locuszoom(assoc_file_name, rangelabel, chromosome,range_start_bp,range_end_bp,pheno_tag, assoc_file_basename)
        loopcount = loopcount + 1
    print('''
*****************************************************************
    ****************************************************
    Termination condition reached, and pc_workhorse
    ended at {0}.
    ****************************************************
******************************************************************
'''.format(time.strftime("%a,%c")))


if __name__ == '__main__':
    main()
    
