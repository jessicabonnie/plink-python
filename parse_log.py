#! /usr/bin/python2.7
#! ./
'''
Created on Oct 20, 2011

@author: Jessica Bonnie
'''
import os
import sys
import getopt
import time
from datetime import datetime, timedelta
from collections import namedtuple
sys.path.insert(0, '/home/jkb4y/h4t1/programs/plink_python/')
import pc_toolbox
import fix_it



global log_folder, summary_name, run_info, freq_loc

#***************************************************
#DEFAULTS - constants used to default the globals

#Defines the path location of the logs folder --->> commandline flag: logfolder
LOG_FOLDER = '/home/jkb4y/work/results/UK_12212011/CC/logs'

#Defines name of summary file --->> commandline flag: summary
SUMMARY_NAME = 'log_summary'

#Defines name of run info file --->> commandline flag: runinfo
RUN_INFO = 'run_info'

#Defines path location of allele frequency file
FREQ_LOC = '/home/jkb4y/work/data/UK_12212011/UK_control.freq'

#Defines path location of original, un map_adapted map file
MAP_LOC = '/home/jkb4y/work/data/UK_12212011/copy_of_UK.bim'

BUILD = 'hg19'

#***************************************************

def three_values(line):
    '''
    Retrieves three values from a specially marked line of the log file.
    Args:
        line -- a specially indicated line from the log file
    Returns:
        last_three -- a list of the three relevant values
    '''
    pieces = line.split()
    last_three = pieces[1:]
    return last_three

def mb_range(line):
    '''Retrieves range information from the log file.
    Args:
        line -- a specially indicated line from the log file
    Returns:
        chr_ref_start_end -- list of four values [chromosome, gene_symbol, start(Mb), end(Mb)]

    '''
    RangeInfo = namedtuple('RangeInfo', 'band,chro,gene,start,end,flag')
    pieces = line.split()
    range_info = RangeInfo._make(pieces[1:])
    print range_info
    #chr_ref_start_end_flag = pieces[1:]
    #return chr_ref_start_end_flag
    return range_info

def snp_p_loop(line):
    '''Retrieves information about the snp, p-value, and loop from log file.
    Args:
        line -- a specially indicated line from the log file
    Returns:
        s_p_l -- list of three values [SNP, p-value, loopcount]

    '''
    pieces = line.split()
    s_p_l = pieces[1:]
    return s_p_l


def one_value(line):
    '''Retrieve a single value from a line of the log file.
    Args:
        line -- a specially indicated line of the log file
    Returns:
        datum -- value from the log file (originally preceeded by indicator)

    '''
    pieces = line.split()
    datum = pieces[1]
    return datum


def read_log(log_summary,log_file,freq_loc,annot_dict):
    '''
    Extract relevant information from log file and append it to a summary file.
    Args:
        log_file -- path to log file which has been formated in a specific manner
        log_summary -- path to which summary information should be written
    Returns:
        Nothing.
    
    '''
    #values used to reinitialize non-range information in output_list
    BLANK = '--'
    empty_list = []
    condition_list = empty_list
    reset_val = [empty_list,BLANK,BLANK,BLANK,BLANK,BLANK,BLANK,BLANK,BLANK,BLANK,'--']
    blank_dict = {'chr':'chr','band':'band','ref':'refgene','start':'start','end':'end','flag':'--',
                  'clist':BLANK,'im':BLANK,'sig':BLANK,'pval':BLANK, 'OR':BLANK,'lo':BLANK,
                  'hi':BLANK,'t':BLANK,'aa1':BLANK,'maf':BLANK,'ma1':BLANK,'ma2':BLANK,
                  'lz':BLANK,'ca':BLANK,'total':'--'}
    order_list = ['chr','band','ref','start','end','flag','clist','im','sig','pval','OR','lo',
                  'hi','t','aa1','maf','ma1','ma2','lz','ca','total']
    reset_list = ['clist','im','sig','pval','OR','lo',
                  'hi','t','aa1','maf','ma1','ma2','lz','ca']
    output_dict = blank_dict
    with open(log_file, mode='r')as log:
        
        with open(log_summary, mode='a')as summary:
            
            for logline in log:
                #if logline begins '$$$', it contains the range info
                if logline.startswith('$$$'):
                    range_info = mb_range(logline)
                    output_dict['band']=range_info.band
                    output_dict['chr']=range_info.chro
                    output_dict['ref']=range_info.gene
                    output_dict['start']=range_info.start
                    output_dict['end']=range_info.end
                    output_dict['flag']=range_info.flag
                if logline.startswith('<<<'):
                    for key in reset_list:
                            output_dict[key]=BLANK
                #if logline begins '%%%', info pertains to loop through plink_association.py
                if logline.startswith('%%%'):
                    #if there is unprinted (and thus unwiped) loop information
                    #in output_dict, LocusZoom did not spit out
                    #a 'Found: ' phrase with reference SNP name. This is bad news!
                    if not output_dict['sig']==BLANK:
                        output_dict['lz'] = 'ERROR'
                        print('ERROR: LZ did not announce a reference SNP!')
                        print output_dict
                        dict_list = []
                        for key in order_list:
                            dict_list.append(output_dict[key])
                        print dict_list
                        summary.write('\t'.join(dict_list)+'\n')
                    #gather info pertaining to this loop through plink
                    snp_pv_loop = snp_p_loop(logline)
##                    print('SNP_PV_LOOP is: ')
##                    print snp_pv_loop
                    output_dict['sig']=snp_pv_loop[0]
                    output_dict['pval']=snp_pv_loop[1]
                    output_dict['ca']=snp_pv_loop[2]
                    #look up the allele frequency of the snp in the freq file
                    maf,ma1,ma2 = pc_toolbox.retrieve_freq(freq_loc, snp_pv_loop[0])
                    output_dict['maf']=maf
                    output_dict['ma1']=ma1
                    output_dict['ma2']=ma2
                    #look up the original Immunochip SNP name in the map file
                    #im = pc_toolbox.retrieve_im(map_loc, snp_pv_loop[0], repair_loc)
                    im = annot_dict[snp_pv_loop[0]].name
                    output_dict['im']=im
                #if logline begins with Odds Ratio etc.
                elif logline.startswith('Regression Coefficient') or logline.startswith('Odds Ratio'):
                    orbeta = logline.strip().partition(':')[2]
                    output_dict['OR']=orbeta
                elif logline.startswith('Coefficient T-Statistic'):
                    stat = logline.strip().partition(':')[2]
                    output_dict['t']=stat
                elif logline.startswith('A1:'):
                    aa1 = logline.strip().partition(':')[2]
                    output_dict['aa1']=aa1
                #if logline begins 'Found: ', obtain name of LZ's reference SNP
                elif logline.startswith('Found:'):
                    lz_snp = one_value(logline)
                    output_dict['lz']=lz_snp
                    if len(condition_list)==0:
                        output_dict['clist']=BLANK
                    else:
                        output_dict['clist']='['+','.join(condition_list)+']'
                    dict_list = []
                    for key in order_list:
                        dict_list.append(output_dict[key])
                    print dict_list
                    summary.write('\t'.join(dict_list)+'\n')
                    condition_list = list()
                    #if output_dict['total'] is '--',
                    #then this is NOT last loop through plink_association.py.
                    #(if it were, p_a wouldn't give a ref SNP, and LZ would,
                    #leading to an expected mismatch.)
                    if output_dict['total']=='--':
                        if not output_dict['lz']==output_dict['sig']:
                            print('ERROR: plink_association and Locus Zoom did not identify the same reference SNP!')
                            summary.write('ERROR: PLINK_ASSOCIATION AND LOCUS ZOOM IDENTIFIED DIFFERENT REFERENCE SNPs! \n')
                    output_dict['total']='--'
                    for key in reset_list:
                            output_dict[key]=BLANK
                    
                #if logline begins ':::', this is final loop through p_a
                elif logline.startswith(':::'):
                    final_loopcount = one_value(logline)
                    output_dict['total']=final_loopcount
                    output_dict['ca']=final_loopcount
                #if logline begins '&&&', it contains a snp in condition list
                elif logline.startswith('&&&'):
                    condition_list.append(one_value(logline))
                    print(condition_list)
                elif logline.startswith('Confidence Interval'):
                    want = logline.partition(':')[2].strip()
                    lo = want.partition('-')[0].strip()
                    hi = want.partition('-')[2].strip()
                    output_dict['lo']=lo
                    output_dict['hi']=hi
                    
                next
                
def summarize_folder(log_folder, freq_loc, summary_name,run_info,annot_dict):
    '''
    Extract summary information from all logfiles in a folder.

    Args:
        folder -- path to logs folder
        summary_name -- desired name of summary file (NO extension)
        freq_loc -- path to allele freqency file
        run_info -- file to contian run info
    Returns:
        Nothing
    
    '''
    summary_file = summary_name + '.tbl'
    run_info_loc = run_info + '.tbl'
    #make a list of all files in the folder with '.log' extension
    log_list = filter(lambda x: x.endswith('.log'),os.listdir(log_folder))
    #if summary file already exists, clear it
    if os.path.exists(summary_file):
        clear = open(summary_file,mode='w')
        clear.close()
    #if snp_count file already exists, clear it
    if os.path.exists(run_info_loc):
        wipe = open(run_info_loc,mode='w')
        wipe.close()
    #create the title rows
    title_list = ['CHR','CHR_BAND','GENE_SYMBOL','START(Mb)','END(Mb)',
                  'FLAG','CONDITION_LIST','IM_CHIP_NAME','SIGNIFICANT_SNP',
                  'P-VALUE','OR/BETA','CI_LO','CI_HI','STAT','LOGISTIC_A1',
                  'MAF','ALLELE_A1','ALLELE_A2','LOCUSZOOM_SNP','COND_ANALYSIS', 'TOTAL_SIG_SNP']
    title_line = '\t'.join(title_list)
    SNP_title_list = ['CHR','CHR_BAND','GENE_SYMBOL','START(Mb)','END(Mb)','#SNPS','BEGIN','END','DURATION']
    SNP_title_line = '  \t'.join(SNP_title_list)
    with open(summary_file, mode='a') as summary:
        #write title row to first line of summary file
        summary.write(title_line+'\n')
    with open(run_info_loc,mode='a') as snippy:
        #write title row of run_info
        snippy.write(SNP_title_line + '\n')
    #append summary of each log file in folder to the summary file
    for log_file in log_list:
        print('Now opening {0} for summarizing.'.format(log_file))
        read_log(summary_file,log_file,freq_loc, annot_dict)
        count_SNPs(log_file,run_info_loc)

        
def usage():
    '''
    Usage information for parse_log.py, and current default information.
    To access:
    python parse_log.py -h  OR  python parse_log.py --help
    '''
    print('''
USAGE: parse_log.py [FLAG] OBJECT
      FLAG              DESCRIPTION                                         CURRENT DEFAULT
    -l, --logfolder     path of folder containing log files                 {0}
    -f, --freq          path of allele freqency file                        {1}
    -s, --summary       desired name of summary file (no extension)         {2}
    -r, --runinfo       desired name of file containing SNP counts (no ext) {3}
    -b, --build         build of rs snps and positions in logs              {4}
    -h, --help          display this usage string
       
    '''.format(LOG_FOLDER,FREQ_LOC,SUMMARY_NAME,RUN_INFO, BUILD))

def cl_arguments(argv):
    '''
    Make use of command line arguments to assign globals.
    Args:
        argv -- sys.argv[1:], sent from main()
    '''
    global summary_name, log_folder, run_info, freq_loc, build
    summary_name = SUMMARY_NAME
    log_folder = LOG_FOLDER
    run_info = RUN_INFO
    freq_loc = FREQ_LOC
    build = BUILD
##    map_loc = MAP_LOC
##    repair_loc = None
    
    try: 
        opts, args = getopt.getopt(argv, "hl:s:r:f:b:",
                                   ["help","logfolder=","summary=",
                                    "runinfo=","freq=","build="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-l","--logfolder"):
            log_folder = arg
        elif opt in ("-f","--freq"):
            freq_loc = arg
        elif opt in ("-b","--build"):
            build = arg
##        elif opt in ("-m","--map"):
##            map_loc = arg
        elif opt in ("-s","--summary"):
            summary_name = arg
        elif opt in ("-r","--runinfo"):
            run_info = arg
##        elif opt in ("--repair-dict"):
##            repair_loc = arg


def count_SNPs(logfile,run_info_loc):
    '''
    Extract information pertaining to the first loop through PLINK,
    e.g. number of SNPs identified in region, duration of association test.
    Purely for use by programmer.

    Args:
        logfile -- location of logfile
    Returns:
        Nothing
    '''
    reset_list = ['--','--','--','--','--','--','--']
    start_time =[]
    snp_num = 0
    end_time = []
    output_list = reset_list
    #open the logfile for reading
    with open(logfile, mode='r') as log:
        #open file to record PLINK run information
        with open(run_info_loc,mode='a'):
            for logline in log:
                #if logline begins with '$$$' it contains range information
                if logline.startswith('$$$'):
                    chr_start_end = mb_range(logline)
                    output_list[:3] = chr_start_end
                #if logline begins 'Analysis started: ' it contains start time
                elif logline.startswith('Analysis started: '):
                    time_tup = logline.partition(':')
                    #time/date are given in asctime, which is not helpful
                    start_time = datetime.strptime(time_tup[2].strip(),"%a %b %d %H:%M:%S %Y")
                    output_list[4] = str(datetime.strftime(start_time,"%H:%M,%m/%d"))
                #if logline begins ... line contains number of SNPs found in region
                elif logline.startswith('After frequency and genotyping pruning, there are '):
                    snp_line_list = logline.split()
                    snp_num = snp_line_list[7]
                    output_list[3]= snp_num
                #if line begins ..., line contains end time
                elif logline.startswith('Analysis finished: '):
                    time_tup = logline.partition(':')
                    #time/date are given in asctime, which is not helpful
                    end_time = datetime.strptime(time_tup[2].strip(),"%a %b %d %H:%M:%S %Y")
                    delta = end_time - start_time
                    output_list[5] = str(datetime.strftime(end_time,"%H:%M,%m/%d"))
                    output_list[6] = str(delta)
                    #once endtime is recorded, first PLINK run is complete, so record info to snippy
                    with open(run_info_loc,mode='a') as snippy:
                        snippy.write(('\t'.join(output_list))+'\n')
                    return

def main():
    
    global log_folder, summary_name, freq_loc, run_info, build #map_loc #, repair_loc
    cl_arguments(sys.argv[1:])
    if log_folder is None:
        usage()
        sys.exit(2)
    os.chdir(log_folder)
    annot_dict_loc = fix_it.locate_annot_dict(build)
    annot_dict = fix_it.build_annot_dict('LOG',annot_dict_loc)
    summarize_folder(log_folder,freq_loc, summary_name,run_info, annot_dict)
    

if __name__=='__main__':
    main() 
