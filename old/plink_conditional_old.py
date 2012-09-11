#! /usr/bin/python2.7
#! ./
'''
@author: Jessica Bonnie
Created on October 3, 2011

'''
import os
import sys
import getopt
import time

global user_script_loc  #location of user's supplied script

global output_folder    #user's supplied output folder
global out_flag         #user's supplied flag for naming files

global p_bound          #user's supplied lowest P-value of interest
global loopagain        #boolean indicating whether another loop will follow current one
global maxloops         #user's supplied maximum number of times to loop through plink
global plink_test       #user's supplied choice of plink association test
global range_start_bp   #translation of user's supplied start to range on chromosome
global range_end_bp     #translation of user's supplied end to range on chromosome
global refgene          #user's supplied reference gene name (used in naming)
global condition_list   #user's supplied location of condition-list

global pheno_loc        #user's supplied location of phenotype file


#*****************************************************
#DEFAULT VALUES - USER MAY ADAPT THESE AS DESIRED

#Defines highest p-value of interest --->> commandline flag: pbound
P_BOUND = 1e-6

#Defines output folder to which to write results --->> commandline flag: outfolder
OUTPUT_FOLDER = None

#Defines maximum number of times to loop thru PLINK --->> commandline flag: loop
MAXLOOPS = 100

#Defines location of user's basic PLINK instruction script --->> commandline flag: script
USER_SCRIPT_LOC = None

#Define PLINK association test to run on data --->> commandline flag: test
PLINK_TEST = 'logistic'
#*******************************************************************

#DEFAULT VALUES - ADAPTATION NOT RECOMMENDED, BUT NOT FORBIDDEN

#Define chromosome on which to run test --->> commandline flag: chromosome
#(if changed - original value to return to: None)
CHROMOSOME = None

#Define flag to add to output file names --->> commandline flag: flag
#(if changed - original value to return to: None)
OUT_FLAG = None

#Define reference gene name to use in naming ---> commandline flag: refgene
#(if changed - original value to return to: None)
REFGENE = None

#Define loc of pre-prepared list of condition SNPS --->> commandline flag: condition-list
#(if changed - original value to return to: None)
CONDITION_LIST = None

#Define loc of phenotype file
PHENO_LOC = None

C_INTERVAL = '0.95'


def log_folder(output_folder):
    log_folder = os.path.join(output_folder, 'logs')
    if not os.path.exists(log_folder):
        os.makedirs(log_folder)
    return log_folder

def convert_bp_to_mb(bp):
    mb = str(round((float(bp)*1e-6),2))
    return mb

def pheno_provided(pheno_loc):
    pheno_list = []
    with open(pheno_loc, mode="r") as phenofile:
        for line1 in phenofile:
            line_split = line1.strip().split()
            pheno_list = line_split[2:]
            return pheno_list

def usage():
    print(

    ('''
USAGE: plink_conditional.py [FLAG] OBJECT
      FLAG                  DESCRIPTION                                     CURRENT DEFAULT
    -o, --outfolder         path of folder where results are to be written  {0}
    -f, --flag              flag to add to output file names                {1}
    -c, --chromosome        number of chromosome on which to run test       {2}
        --from-bp	    start of range on chromosome in bp		    
	--to-bp		    end of range on chromosome in bp		    
	--from-kb           start of range on chromosome in kb
        --to-kb             end of range on chromosome in kb
        --from-mb           start of range on chromosome in mb
        --to-mb             end of range on chromosome in mb
	--refgene	    Gene Name of Interest (for naming purposes)	    {3}
    -p, --pbound            highest acceptable p-value (cut-off point)      {4}
    -l, --loop              maximum number of loops through PLINK           {5}
    -t, --test              which association test PLINK should run         {6}
    -s, --script            pathname of script to feed into PLINK           {7}
        --condition-list    pathname of previously compiled condition-list
        --pheno             pathname of plink pheno file
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
    ).format(OUTPUT_FOLDER,OUT_FLAG,CHROMOSOME,REFGENE,
             str(P_BOUND),str(MAXLOOPS),PLINK_TEST,USER_SCRIPT_LOC))

def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global output_folder, plink_test, p_bound, chromosome
    global maxloops, out_flag, user_script_loc, pheno_loc
    global range_start_bp, range_end_bp, refgene, condition_list,c_interval
    
    
    plink_test = PLINK_TEST
    p_bound = float(P_BOUND)
    chromosome = CHROMOSOME
    output_folder = OUTPUT_FOLDER
    out_flag = OUT_FLAG
    user_script_loc = USER_SCRIPT_LOC
    maxloops = MAXLOOPS
    refgene = REFGENE
    condition_list = CONDITION_LIST
    pheno_loc = PHENO_LOC
    c_interval = None
    
    try: 
        opts, args = getopt.getopt(argv, "ho:s:f:t:c:p:l:g:",
                                   ["help","outfolder=","script=","flag=","test=",
                                    "chromosome=","pbound=","loop=","refgene=",
                                    "from-kb=","to-kb=","from-mb=","to-mb=",
				    "from-bp=","to-bp=","condition-list=",
                                    "pheno=","ci="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-f","--flag"):
            out_flag = arg
        elif opt in ("-s","--script"):
            user_script_loc = arg
        elif opt in ("-o","--outfolder"):
            output_folder = arg
        elif opt in ("-t","--test"):
            plink_test = arg
        elif opt in ("-c","--chromosome"):
            chromosome = arg
        elif opt in ("--from-kb"):
            range_start_bp = int(1000 * float(arg))
        elif opt in ("--to-kb"):
            range_end_bp = int(1000 * float(arg))
        elif opt in ("--from-mb"):
            range_start_bp = int(1e6 * float(arg))
        elif opt in ("--to-mb"):
            range_end_bp = int(1e6 * float(arg))
	elif opt in ("--from-bp"):
	    range_start_bp = int(arg)
	elif opt in ("--to-bp"):
	    range_end_bp = int(arg)
        elif opt in ("-g","--refgene"):
            refgene = arg
        elif opt in ("-p","--pbound"):
            p_bound = float(arg)
        elif opt in ("-l","--loop"):
            maxloops = int(arg)
        elif opt in ("--condition-list"):
            condition_list = arg
        elif opt in ("--pheno"):
            pheno_loc = arg
        elif opt in ("--ci"):
            c_interval = arg
            
def create_command(script_loc , pheno_tag=None):
    global chromosome,range_start_bp,range_end_bp,plink_test
    global p_bound, pheno_loc, c_interval
    global output_folder, maxloops, out_flag, condition_list
    cmd_list = ["python pc_workhorse.py","--chromosome",chromosome,
                "--from-bp", str(range_start_bp),"--to-bp", str(range_end_bp),
                "--pbound",str(p_bound), "--test", plink_test,
                "--script", script_loc, "--outfolder",output_folder,
                "--loop",str(maxloops)]
    of = ""
    if out_flag is not None:
        of = '_'+out_flag
        f_cmd = ["--flag",out_flag]
        cmd_list.extend(f_cmd)
        
    logname = ""
    
    log_start = 'Chr'+chromosome+'_'
    if refgene is None:
        r_start = convert_bp_to_mb(range_start_bp)
        r_end = convert_bp_to_mb(range_end_bp)
        logname =log_start + r_start + '-'+r_end+of
    else:
        r_cmd = ["--refgene",refgene]
        cmd_list.extend(r_cmd)
	logname = log_start + refgene+of
    
    lf = log_folder(output_folder)
    #log_file = os.path.join(lf,log)
    if pheno_tag is not None:
        logname = logname + '_'+pheno_tag
        p_cmd = ["--pheno-tag",pheno_tag]
        cmd_list.extend(p_cmd)
    ci_cmd = []
    if c_interval is None and plink_test in ('logistic',
                                             'linear','assoc',
                                             'fisher'):
        c_interval = C_INTERVAL
    if c_interval is not None:
        ci_cmd = ["--ci",str(c_interval)]
        cmd_list.extend(ci_cmd)
    
    t_cmd = []
    
    if condition_list is None:
        log = logname + '.log'
        log_file = os.path.join(lf,log)
        #t_cmd = ["2>&1|tee",log_file]
    else:
        c_cmd = ["--condition-list",condition_list]
        cmd_list.extend(c_cmd)
        log = logname + '~.log'
        log_file = os.path.join(lf,log)
	#t_cmd = ["2>&1|tee -a",log_file]

    t_cmd = ["2>&1|tee",log_file]
    cmd_list.extend(t_cmd)
    
    cmd = ' '.join(cmd_list)
    print ('My log file is {0} and my command is {1}'.format(log_file,cmd))
    return cmd

def main():
    global chromosome, range_start_bp, range_end_bp, pheno_loc, user_script_loc
    
    cl_arguments(sys.argv[1:])
    if chromosome is None:
        print('''
    BAD NEWS BEARS! You have not specified a chromosome on which to run
    your conditional analysis. Please try again!
              ''')
        sys.exit(2)
    chr_range = 'chr'+ chromosome+':'+str(int(range_start_bp))+'-'+str(int(range_end_bp))

    if pheno_loc is not None:
        pheno_list = pheno_provided(pheno_loc)
        base, ext = os.path.splitext(user_script_loc)
        
        for pheno_tag in pheno_list:
            pheno_script = base + '_'+pheno_tag+ext
            with open(pheno_script, mode="w") as phephe:
                phephe.write("--pheno " + pheno_loc + '\n')
                phephe.write("--pheno-name "+pheno_tag + '\n')
                with open(user_script_loc, mode="r") as skippy:
                    for skip in skippy:
                        phephe.write(skip)
            user_script = pheno_script
            cmd = create_command(user_script, pheno_tag=pheno_tag)
            os.system(cmd)
    else:
        cmd = create_command(user_script_loc)
        os.system(cmd)
    
if __name__ == '__main__':
    main()
