#! /usr/bin/python2.6
#! ./
'''
Created on Jun 26, 2012

@author: Jessica Bonnie
'''
import os
import subprocess
import sys
import getopt
from collections import namedtuple
import shutil

import fix_it
import pc_toolbox
import testsigs_wrap


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
            uss = uss + line
    return uss


def write_script(outfile, new_script_loc, user_script_loc, condition_snp, single, snpstar):
    global chromosome, range_start_bp
    global range_end_bp
    user_script_string = get_user_script_text(user_script_loc)
    with open(new_script_loc, mode="w") as new_script:
        new_script.write(user_script_string + '\n')
        new_script.write('--chr ' + chromosome +'\n'+
                         '--from-bp ' + str(range_start_bp) +'\n'+
                         '--to-bp ' + str(range_end_bp) +'\n'+
                         '--out ' + outfile+'\n'+
                         '--condition '+condition_snp)
        if single:
            new_script.write('\n--snps '+snpstar+','+condition_snp)

        

def make_snp_pos_list(input_assoc, output_name):
    snplist = list()
    line1 = True
    with open(input_assoc, mode="r") as assoc:
        for aline in assoc:
            if line1:
                index_dict = pc_toolbox.read_assoc_titles(aline,
                                                          .95, 'logistic')
                line1 = False
            else:
                lsplit = aline.strip().split()
                snp_tuple = (lsplit[index_dict['snp']], lsplit[index_dict['pos']])
                snplist.append(snp_tuple)
    with open(output_name, mode="w") as output:
        for snp in snplist:
            output.write('\t'.join(snp)+'\n')
    return snplist

def filter_result(result_file, snpstar):
    line1 = True
    with open(result_file, mode="r") as result:
        for line in result:
            if line1:
                index_dict = pc_toolbox.read_assoc_titles(line,
                                                          .95, 'logistic')
                line1 = False
            else:
                lsplit = line.strip().split()
                if lsplit[index_dict['snp']] == snpstar:
                    #print('$$$ '+ line.strip())
                    info_tuple = pc_toolbox.assoc_tuple(line, index_dict)
                    #print info_tuple
                    return info_tuple

def summarize_table(table_loc, summary_loc, snpstar):
    line1 = True
    pbound = .05
    rbound = .8
    PR, Pr, pR, pr = 0, 0, 0, 0
    NAr = 0
    NAR = 0
    NAR5 = 0
    NAq = 0
    pq = 0
    Pq = 0
    PR5 = 0
    pR5 = 0
    snpi = None
    with open(table_loc, mode="r")as table:
        for line in table:
            lsplit = line.strip().split()
            if line1:
                r2_i = lsplit.index('r2')
                p_i = lsplit.index('SNP*_pvalue')
                snpi_i = lsplit.index('conditional_snp')
                line1 = False
            elif not snpstar == lsplit[snpi_i]:
                if lsplit[p_i] == 'NA':
                    if lsplit[r2_i] == '???':
                        NAq = NAq + 1
                    else:
                        if float(lsplit[r2_i]) > rbound:
                            NAR = NAR + 1
                        elif float(lsplit[r2_i]) > .5:
                            NAR5 = NAR5 + 1
                        else:
                            NAr = NAr + 1
                elif float(lsplit[p_i]) > pbound:
                    if lsplit[r2_i] == '???':
                        Pq = Pq + 1
                    else:
                        if float(lsplit[r2_i]) > rbound:
                            PR = PR + 1
                        elif float(lsplit[r2_i]) > .5:
                            PR5 = PR5 + 1
                        else:
                            Pr = Pr + 1
                else:
                    if lsplit[r2_i] == '???':
                        pq = pq + 1
                    else:
                        if float(lsplit[r2_i]) > rbound:
                            pR = pR + 1
                        elif float(lsplit[r2_i]) > .5:
                            pR5 = pR5 + 1
                        else:
                            pr = pr + 1
##
##                else:
##                    r2 = float(lsplit[r2_i])
##                if lsplit[p_i] == 'NA':
##                    p = pbound + .01
##                else:
##                    p = float(lsplit[p_i])
##                if p > pbound:
##                    if r2 > rbound:
##                        gpgr = gpgr + 1
##                    else:
##                        gplr = gplr + 1
##                else:
##                    if r2 > rbound:
##                        lpgr = lpgr + 1
##                    else:
##                        lplr = lplr + 1
    with open(summary_loc, mode = "w") as summary:
        summary.write('\t'.join(['.','p>0.05','p<=0.05','p=NA'])+'\n')
        summary.write('\t'.join(['r2>0.8',str(PR),str(pR),str(NAR)])+'\n')
        summary.write('\t'.join(['0.5<r2<=0.8',str(Pr),str(pr),str(NAr)])+'\n')
        summary.write('\t'.join(['r2<=0.5',str(PR5),str(pR5),str(NAR5)])+'\n')
        summary.write('\t'.join(['r2=???',str(Pq),str(pq),str(NAq)])+'\n')

        

def usage():
    print("usage goes here")

def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global outfolder, assoc, chromosome,snpstar, ldfolder, freq_loc
    global out_flag, user_script_loc, hitstring
    global range_start_bp, range_end_bp, hit1
    global region_id, build, single, multi, hit_index
    
    out_flag = ''
    plink_test = None
    region_id = None
    build = 'hg19'
    assoc = None
    snpstar = None
    outfolder = None
    range_start_bp = None
    range_end_bp = None
    user_script_loc = None
    single = False
    ldfolder = testsigs_wrap.LDFOLDER
    freq_loc = testsigs_wrap.FREQ_LOC
    multi = False
    hit_index = None
    hit1 = False
    hitstring = None
    
    try: 
        opts, args = getopt.getopt(argv, "ho:s:f:c:",
                                   ["help","outfolder=","script=","flag=",
                                    "chromosome=",
                                    "start=","end=",
                                    "assoc=","snpstar=","single",
                                    "region-id=","build=","interrupt",
                                    "freq=","ldfolder=","multi=","hit1="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-f","--flag"):
            out_flag = '_'+arg
        elif opt in ("-s","--script"):
            user_script_loc = arg
        elif opt in ("-o","--outfolder"):
            outfolder = arg
            print("Read folder arg!")
        elif opt in ("-c","--chromosome"):
            chromosome = arg
        elif opt in ("--snpstar"):
            snpstar = arg
        elif opt in ("--ldfolder"):
            ldfolder = arg
        elif opt in ("--freq"):
            freq_loc = arg
        elif opt in ("--single"):
            single = True
        elif opt in ("--multi"):
            multi = True
            hit_index = arg
        elif opt in ("--hit1"):
            hit1 = True
            hitstring = arg
        elif opt in ("--assoc"):
            assoc = arg
        elif opt in ("--start"):
            range_start_bp = int(1e6 * float(arg))
        elif opt in ("--end"):
            range_end_bp = int(1e6 * float(arg))
        elif opt in ("--region-id"):
            region_id = arg
        elif opt in ("--build"):
            build = arg
            print("BUILD IS SET TO: {0}!".format(arg))



def main(argv):
    global outfolder, assoc, chromosome,snpstar
    global out_flag, user_script_loc
    global range_start_bp, range_end_bp, hit_index, hitstring
    global region_id, build, single, ldfolder, freq_loc, multi, hit1
    
    cl_arguments(argv)
    placeholder = ''
    annot_dict_loc = fix_it.locate_annot_dict(build)
    annot_dict = fix_it.build_annot_dict('LOG',annot_dict_loc)
    snpstar_im = annot_dict[snpstar].name
    table_folder = os.path.join(outfolder, 'ResultTables')
    summary_folder = os.path.join(outfolder, 'SummaryTables')
    str_hit_index = str(hit_index)
    if hit_index < 10:
        str_hit_index = '0'+str_hit_index
    if not os.path.exists(table_folder):
        os.makedirs(table_folder)
    if not os.path.exists(summary_folder):
        os.makedirs(summary_folder)
    chr_folder = os.path.join(outfolder, 'chr{0}'.format(chromosome))
    if not os.path.exists(chr_folder):
        os.makedirs(chr_folder)
    reg_folder = os.path.join(chr_folder, region_id)
    if multi:
        placeholder = '_'+snpstar+'_'+str_hit_index
    if not os.path.exists(reg_folder):
        os.makedirs(reg_folder)
    super_outbase = os.path.join(reg_folder, region_id)
    list_loc = super_outbase + '_snps.list'
    
    table_loc = os.path.join(chr_folder, region_id+placeholder+'.tbl')
    new_table_loc = os.path.join(table_folder, region_id+placeholder+'.tbl')
    summary_table_loc = os.path.join(summary_folder, region_id+placeholder+'.tbl')
    print(table_loc)
    snplist = make_snp_pos_list(assoc,list_loc)
    ld_loc = os.path.join(ldfolder,'chr{0}'.format(chromosome),
                              '{0}_r2_0.ld'.format(region_id))
    table = open(table_loc, mode="w")
    table.write('\t'.join(['SNP*','SNP*_pos','SNP*_im','conditional_snp',
                           'csnp_im','csnp_pos','SNP*_pvalue','OR','ci_lo','ci_hi',
                           'a1',"r2","csnp_freq","csnp_freq_a1"])+'\n')
    table.close()
    index = 1
    for snp_tuple in snplist:
        snp = snp_tuple[0]
        snp_pos = snp_tuple[1]
        snp_im = annot_dict[snp].name
        corrected_snp = snp.replace(':','_')
        outbase = super_outbase+'_'+corrected_snp
        script_loc = outbase + '.script'
        assoc_out = outbase +'.assoc.logistic'
        if hit1:
            write_script(outbase,script_loc, user_script_loc, snp, single, hitstring)
            plink(script_loc)
            print('''
**********************************************************************
    ************************************************************
    The data will be conditioned on the following SNP:
    %%%             {0}
    This is snp #{1} in a list of {2}.
        ******************************************************
**********************************************************************
'''.format(snp,index,len(snplist)))
##            if snp == snpstar:
##                write_script(outbase,script_loc, user_script_loc, snp, single, hitstring)
##                plink(script_loc)
        elif not multi:
            write_script(outbase,script_loc, user_script_loc, snp, single, snpstar)
            plink(script_loc)
            print('''
**********************************************************************
    ************************************************************
    The data will be conditioned on the following SNP:
    %%%             {0}
    This is snp #{1} in a list of {2}.
        ******************************************************
**********************************************************************
'''.format(snp,index,len(snplist)))
            
        info = filter_result(assoc_out,snpstar)
        snp_freq = pc_toolbox.retrieve_freq(freq_loc, snp)
        r2 = pc_toolbox.retrieve_r2(snpstar,snp,ld_loc)
        if r2 is None:
            r2 = "???"
        
        index = index + 1
        with open(table_loc, mode='a') as table:
            table.write('\t'.join([snpstar,info.pos, snpstar_im,
                                    snp,snp_im,snp_pos,str(info.p),
                                    info.OR,info.lo,info.hi,
                                    info.a1,r2,snp_freq[0],snp_freq[1]])+'\n')
    shutil.copy(table_loc, new_table_loc)
    summarize_table(new_table_loc, summary_table_loc, snpstar)




if __name__ == '__main__':
    main(sys.argv[1:])
