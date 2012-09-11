#! /usr/bin/python2.7
'''
Created on Sep 14, 2011

@author: Jessica Bonnie
'''
import os
import sys
import getopt
import subprocess
from collections import namedtuple
import fix_it

OUTFOLDER = '/home/jkb4y/work/results/HapMap/DupConc/'
LOG_LOC = '/home/jkb4y/work/results/HapMap/03162012-Complete_HapMap_Immuno_Full_Data_Table.log'
BFILE = '/home/jkb4y/work/results/HapMap/03162012-Complete_HapMap_Immuno_Full_Data_Table'

def plink_recode(bfile, recode_loc, out, script_loc):
    with open(script_loc, mode="w") as ps:
        ps.write("--bfile "+bfile+"\n")
        ps.write("--update-ids "+recode_loc+"\n")
        ps.write("--out "+ out + "\n")
        ps.write("--make-bed"+"\n")
        ps.write("--noweb")

def plink_remove(bfile, remove_loc, out, script_loc):
    with open(script_loc, mode="w") as ps:
        ps.write("--bfile "+bfile+"\n")
        ps.write("--remove "+remove_loc+"\n")
        ps.write("--out "+ out + "\n")
        ps.write("--make-bed"+"\n")
        ps.write("--noweb")

def plink_keep(bfile, keep_loc, out, script_loc):
    with open(script_loc, mode="w") as ps:
        ps.write("--bfile "+bfile+"\n")
        ps.write("--keep "+keep_loc+"\n")
        ps.write("--out "+ out + "\n")
        ps.write("--make-bed\n")
        ps.write("--noweb")

def plink_merge(remove_b, keep_b, out, script_loc):
    with open(script_loc, mode="w") as ps:
        ps.write("--bfile "+remove_b+"\n")
        ps.write("--bmerge "+keep_b+'.bed '+keep_b+'.bim '+keep_b +'.fam\n')
        ps.write("--out "+ out + "\n")
        #ps.write("--make-bed\n")
        ps.write("--noweb\n")
        ps.write("--merge-mode 6")

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


def make_recode(log_loc, outfolder):
    recode1_loc = os.path.join(outfolder,"recode1.txt")
    recode1 = open(recode1_loc, mode = "w")
    recode2_loc = os.path.join(outfolder,"recode2.txt")
    recode2 = open(recode2_loc, mode = "w")
    triple_loc = os.path.join(outfolder,"removed_not_returned.txt")
    trip_file = open(triple_loc, mode = "w")
    trip_list = list()
    with open(log_loc, mode="r") as log:
        for line in log:
            if line.startswith("Duplicate individual found:"):
                print line
                (pre, sep, info) = line.partition(':')
                print (pre, sep, info)
                
                id_pair = tuple(info.strip().strip('[]').split())
                print id_pair
                id_tup = tuple([id_pair[0],id_pair[1] + '.dup'])
                (fam, new_samp) = id_tup
                if id_tup in trip_list:
                    trip_file.write('\t'.join([id_pair[0],id_pair[1]])+'\n')
                    continue
##                if tuple(id_tup) in trip_list:
##                    (fam, new_samp) = rename_as_necessary(id_tup, trip_list)
##                else:
##                    (fam, new_samp) = id_tup
                print fam, new_samp
                trip_list.append((fam, new_samp))
                print trip_list
                recode1.write('\t'.join([id_pair[0],id_pair[1],fam, new_samp])+'\n')
                recode2.write('\t'.join([fam, new_samp,id_pair[0],id_pair[1]])+'\n')
        recode1.close()
        recode2.close()
    return recode1_loc, recode2_loc, triple_loc
        

def rename_as_necessary(new_pair,trip_list):
    '''
    Determines if the desired new name for the input map file already exists,
    and, if so, chooses another name.
    Args:
        new_orig -- target to which to map file could be moved
        ext -- extension of the map file
    Returns:
        next_rename -- target name to which the original map file can
                        be saved without overwriting any other file
    '''
    next_pair = new_pair
    print next_pair
    new_samp = new_pair[1]
    if tuple(new_pair) in trip_list:
        print("There's already at least one duplicate of this!")
        if new_samp[-1].isdigit():
            pieces = new_samp.rpartition('.dup')
            up1 = int(pieces[2])+ 1
            next_samp = pieces[0]+'.dup'+ str(up1)
            next_pair = (new_pair[0],next_samp)
            print(next_pair)
            next_pair = rename_as_necessary(next_pair, trip_list)
        else:
            print("I have entered where I should add a 2")
            next_pair = (new_pair[0],new_samp + '2')
            print next_pair
            next_pair = rename_as_necessary(next_pair, trip_list)
    print("I have determined the new pair to be:")
    print next_pair
    return next_pair
                  
                
def summarize_diff(merge_out):
    old_samp = 'IID'
    old_fam = 'FID'
    counter = 1
    line1 = True
    DiffInfo = namedtuple("DiffInfo","count,fam,samp")
    diff_list = list()
    #snp_dict = {}
    diff_file = [line.strip() for line in open(merge_out+'.diff', mode="r")]
    diff_file.sort(key= lambda line: line.split()[2])
    for line in diff_file:
        d = line.strip().split()
        new_samp = d[2]
        new_fam = d[1]
        if not new_samp == 'IID':
            if line1:
                old_samp = new_samp
                old_fam = new_fam
                line1 = False
            elif new_samp == old_samp and new_fam == old_fam:
                counter = counter + 1
            else:
                diff_info = DiffInfo(count=counter,fam=old_fam, samp=old_samp)
                diff_list.append(diff_info)
                    #snp_dict[old_snp]= str(counter)
                old_samp=new_samp
                old_fam=new_fam
                counter = 1
##    with open(merge_out+'.diff', mode="r") as diff:
##        
##        for line in diff:
##            d = line.strip().split()
##            new_samp = d[2]
##            new_fam = d[1]
##            if not new_samp == 'IID':
##                if line1:
##                    old_samp = new_samp
##                    line1 = False
##                elif new_samp == old_samp and new_fam == old_fam:
##                    counter = counter + 1
##                else:
##                    diff_info = DiffInfo(count=counter,
##                                         fam=old_fam, samp=old_samp)
##                    diff_list.append(diff_info)
##                    #snp_dict[old_snp]= str(counter)
##                    old_samp=new_samp
##                    old_fam=new_fam
##                    counter = 1
##        #print snp_dict
##    sort_list = sorted(diff_list, reverse=True)
##    print sort_list
    with open(merge_out+"_DiffSummary_samples.txt",mode="w") as summary:
        summary.write('\t'.join(['FID','IID','COUNT'])+'\n')
        for diff in diff_list:
            summary.write(diff.fam+'\t'+diff.samp + '\t'+str(diff.count)+'\n')                

def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global bfile, outfolder, log_loc
    bfile = BFILE
    outfolder = OUTFOLDER
    log_loc = LOG_LOC
    try: 
        opts, args = getopt.getopt(argv, "hb:o:l:",
                                   ["help","bfile=","outfolder=","log="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-b","--bfile"):
            bfile = arg
        elif opt in ("-l","--log"):
            log_loc = arg
        elif opt in ("-o","--outfolder"):
            outfolder = arg
        
def usage():
    print('''
USAGE: dup_check.py [FLAG] OBJECT
      FLAG                  DESCRIPTION                         CURRENT DEFAULT
    -b, --bfile             filepath of plink binary file
    -o, --outfolder         filepath to output folder
    -l, --log               filepath to plink log complaining of duplicates
    -h, --help              display this usage string

''')



def main(argv):
    global bfile, outfolder, log_loc
    cl_arguments(argv)
    out_base = os.path.join(outfolder,'hapmap')
    recode1_loc, recode2_loc, triple_loc = make_recode(log_loc, outfolder)
    keep_bfile = out_base + "_keepDup"
    remove_bfile = out_base + "_removeDup"
    #remove the triples
    bfile_stripped = out_base + '_stripped'
    strip_script = out_base + '_stripScript.txt'
    plink_remove(bfile, triple_loc, bfile_stripped, strip_script)
    plink(strip_script)
    #remove the duplicates:
    remove_script = out_base + '_removeScript.txt'
    plink_remove(bfile_stripped, recode1_loc, remove_bfile,remove_script)
    plink(remove_script)
    
    bim = remove_bfile +'.bim'
    fix_args = ['--map',bim,'--hapmap','--build', 'hg18']
    fix_it.main(fix_args)
    #keep the duplicates
    keep_script = out_base + '_keepScript.txt'
    plink_keep(bfile_stripped, recode1_loc, keep_bfile,keep_script)
    plink(keep_script)

    bim = keep_bfile +'.bim'
    fix_args = ['--map',bim,'--hapmap','--build', 'hg18']
    fix_it.main(fix_args)

     #run concordance check
    merge_script = out_base + '_mergeScript.txt'
    merge_out = out_base + '_DupMerge'
    plink_merge(remove_bfile, keep_bfile,merge_out, merge_script)
    plink(merge_script)

    #summarize the results of the diff file
    summarize_diff(merge_out)
    #fix the map file
##    bim = remove_bfile +'.bim'
##    fix_args = ['--map',bim,'--hapmap','--build', 'hg18']
##    fix_it.main(fix_args)

if __name__=='__main__':
    main(sys.argv[1:])
    

