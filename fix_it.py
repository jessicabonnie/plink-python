#! /ust/bin/python2.7
#! ./
'''
Created Jan 25, 2012

@author: Jessica Bonnie
'''
from collections import namedtuple
import os
import sys
import getopt
import shutil

import meta_toolbox
import pc_toolbox


AnnotInfo = namedtuple("AnnotInfo","name,lz,rs,hg19_chr,hg19_pos,hg18_chr,hg18_pos,band,annot,meta_key")#,meta_key")
ANNOT_LOC = '/home/jkb4y/work/data/quinlan-immunochip-snps-annotated-2012-Aug-03.txt'
#ANNOT_DICT = '/home/jkb4y/work/data/annot.dict'
#THIS SNP needs to be itself in both hg18 and hg19 to avoid duplication
RS_KEEP = ['rs10127859']
#THIS SNP needs to be itself in hg18, but chr:pos in hg19
RS_SWITCH = ['rs10782171']
HG19_CHR_TITLE = 'hg19_chrom'
HG18_CHR_TITLE = 'hg18_chrom'
NAME_TITLE = 'name'
RS_ID_TITLE = 'rsID'
HG18_POS_TITLE = 'hg18_end'
BAND_TITLE = 'band'
HG19_POS_TITLE = 'hg19_end'
ANNOT_TITLE = 'functionGVS'
BUILD = 'hg18'
KEEP_IM_LOC = '/home/jkb4y/work/data/keep_im_list.txt'
DUPLICATE_FIX_DICT = {'chr6_106659585':'rs17066588',
                      'chr2_230785920':'rs11556887',
                      'chr2_230857352':'rs4972946'}
#DUPLICATE_FIX_DICT = {'chr6_106659585':'rs17066588','imm_2_230785920':'rs11556887','imm_2_230857352':'rs4972946'}
IM_DUPLICATES = '/home/jkb4y/work/data/HapMap/hapmap_removeDup~1_imDups.list'
FAM_DIFF = {'imm_2_230857352':'rs4972946'}
#THESE SNPs needs to be chr:pos in hg19, because locuszoom can't "see" it
LZ_CHRPOS = ['imm_20_1502639','imm_12_110591701']

def read_annot(annot_indices, line, im_keep, build,im_duplicates):
    alist = line.strip().split()
    hg18_chr=alist[annot_indices['hg18_chr']].replace('X','23').replace('Y','24').replace('XY','23')

    hg19_chr=alist[annot_indices['hg19_chr']].replace('X','23').replace('Y','24').replace('XY','23')

    rs=alist[annot_indices['rs']]
    hg18_pos=alist[annot_indices['hg18_pos']]
    hg19_pos=alist[annot_indices['hg19_pos']]
    band=alist[annot_indices['band']]
    name=alist[annot_indices['name']]
    annot=alist[annot_indices['annot']]

    #print im_dups
    meta_key = name
    normally = ['imm_2_230785920','imm_2_230857352','imm_6_106659585']
    weimin = ['imm_2_230857352','chr6_106659585']
    quinlan_repair = ['seq-t1d-12-56449572-T-C']
    if name in weimin:
        meta_key = rs
    elif name in RS_KEEP or name in RS_SWITCH:
        meta_key = name
    elif name in im_duplicates:
        meta_key = name
    else:
        if (not rs == '0') and (name not in quinlan_repair):
            meta_key = rs
    #ref = rs

    if build == 'hg19':
        if rs == '0' or name in im_duplicates:# or name in LZ_CHRPOS:
        #if rs == '0' or name in im_keep:
#These SNPs no longer have rs-ids, but cannot become chr:pos because they have a dup
            if name in RS_KEEP:
                lz=name
            elif name in DUPLICATE_FIX_DICT.keys():
                lz = DUPLICATE_FIX_DICT[name]
            else:
                lz = "chr"+hg19_chr+":"+hg19_pos
##        elif name in DUPLICATE_FIX_DICT:
##            lz = DUPLICATE_FIX_DICT[name]
#These SNPs have new RS-ids but cannot use them because their duplicates are
        elif name in RS_SWITCH:
            lz="chr"+hg19_chr+":"+hg19_pos
        else:
            lz = rs
    if build == 'hg18':
        if name.startswith('rs'):
            lz = name
        elif name in DUPLICATE_FIX_DICT.keys():
            lz = DUPLICATE_FIX_DICT[name]
        else:
            lz = "chr"+hg18_chr+":"+hg18_pos
    annot_info = AnnotInfo(hg18_chr=hg18_chr,hg18_pos=hg18_pos,
                           hg19_chr=hg19_chr,hg19_pos=hg19_pos,
                           name=name,rs=rs,
                           lz=lz,band=band,annot=annot,meta_key=meta_key)
##                if not ref in keep_rs_list:
    return annot_info

def read_list_file(list_loc):
    with open(list_loc,mode="r") as listy:
        list_info = list()
        for info in listy:
            i = info.strip()
            list_info.append(i)
    return list_info

def write_annot_dict(annot_loc, build):
    counter = 0
    dict_loc = locate_annot_dict(build, annot_loc)
    annot_dict = open(dict_loc, mode="w")
    im_duplicates = read_list_file(IM_DUPLICATES)
    im_keep = list()
    with open(KEEP_IM_LOC, mode="r") as im:
        for im_line in im:
            im_keep.append(im_line.strip())
    with open(annot_loc, mode = "r") as annot:
        line1 = True
        for aline in annot:
            if line1:
                a_indices = read_annot_titles(aline)
                line1 = False
                annot_dict.write('\t'.join(['MitoA3721G','rs41355750','rs41355750',
                                            'M','3720',
                                            'M','3721',
                                            'band','no-annotation','MitoA3721G']) +'\n')
                print("DID IT")
            else:
                a_info = read_annot(a_indices,aline,im_keep,build,im_duplicates)
                annot_dict.write('\t'.join(a_info) +'\n')
    annot_dict.close()

def locate_fixed_meta(table_loc, build):
    (basepath, ext) = os.path.splitext(table_loc)
    output_loc = basepath + '_lz'+'_'+build+ext
    return output_loc

def get_lz_annot(annot_dict,snp):
    if annot_dict[snp].rs == 'rs6043409':
        return 'nonsyn'
    #there are 2 snps in chr16 marked as nonsyn which we want to draw as non-annot:
    if annot_dict[snp].rs =='rs4381620' or annot_dict[snp].rs =='rs7193871':
        return 'no-annotation'
    annot = annot_dict[snp].annot
    utr_list = ['utr-3','utr-5']
    coding_list = ['coding-synonymous','coding-notMod3','coding-synonymous-near-splice']
    splice_list = ['splice-5','splice-3']
    nonsyn_list = ['missense','stop-gained','stop-lost','missense-near-splice','stop-gained-near-splice']
    if annot in utr_list:
        return 'utr'
    elif annot in splice_list:
        return 'splice'
    elif annot in coding_list:
        return 'coding'
    elif annot in nonsyn_list:
        return 'nonsyn'
    else:
        return 'no-annotation'

def fix_table(table_loc, annot_dict, purpose, build):
    line1 = True
    counter = 0
    if purpose in ['MAP','HAPMAP','FAMMAP']:
##        ##    (basepath, ext) = os.path.splitext(table_loc)
##    if not meta:
        (basepath, ext) = os.path.splitext(table_loc)
        orig_rename = basepath + '~'
        new_loc = str(rename_as_necessary(orig_rename, ext)) + ext
        shutil.copy(table_loc, new_loc)
        output_loc = table_loc
        table_loc = new_loc
    else:
        output_loc = locate_fixed_meta(table_loc, build)
##    else:
##        output_loc = locate_fixed_meta(table_loc, build)
    with open(table_loc, mode="r") as table:
        output = open(output_loc, mode="w")
        index_dict = {'chr':0,'pos':3,'snp':1}
        error_list = list()
        for line in table:
            if counter > 5000:
                print line_list
            if line1 and purpose not in ['MAP','HAPMAP','FAMMAP']:
                if purpose == "META":
                    index_dict = pc_toolbox.read_meta_titles(line)
                if purpose == "FAMILY":
                    index_dict = pc_toolbox.read_fam_titles(line)
##                if meta:
##                index_dict = meta_toolbox.read_meta_titles(line)
                print index_dict
                line_list = line.strip().split()
                if not purpose in ['MAP','HAPMAP','FAMMAP']:
                    line_list.append('annotation')
                output.write('\t'.join(line_list)+'\n')
                index_dict['annot']=len(line_list)-1
                print index_dict
                line1 = False
##                else:
##                    index_dict = {'chr':0,'pos':3,'snp':1}
##                output.write(line)
##                line1 = False
            else:
                line_list = line.strip().split()
                snp = line_list[index_dict['snp']]
                if purpose in ['FAMILY','FAMMAP'] and snp in FAM_DIFF.keys():
                    print("MADE IT INTO CASE")
                    print snp
                    snp = FAM_DIFF[snp]
                    print snp
                if not purpose in ['MAP','HAPMAP','FAMMAP']:
                    lz_annot = get_lz_annot(annot_dict, snp)
                    line_list.append(lz_annot)
##                chro = line_list[index_dict['chr']]
##                pos = line_list[index_dict['pos']]
                if build == 'hg18':
                    try:
                        line_list[index_dict['chr']]= annot_dict[snp].hg18_chr.replace('X','23').replace('Y','24').replace('XY','23').replace('M','26')
                        line_list[index_dict['pos']]= annot_dict[snp].hg18_pos
                    except KeyError:
                        error_list.append(list(line_list))
##                        line_list[index_dict['chr']] = 'ERROR'
##                        line_list[index_dict['pos']]= 'ERROR'
                        
                if build == 'hg19':
                    try:
                        line_list[index_dict['chr']]= annot_dict[snp].hg19_chr.replace('X','23').replace('Y','24').replace('XY','23').replace('M','26')
                        line_list[index_dict['pos']]= annot_dict[snp].hg19_pos
                    except KeyError:
                        error_list.append(list(line_list))
##                        line_list[index_dict['chr']] = 'ERROR'
##                        line_list[index_dict['pos']]= 'ERROR'
##                line_list[index_dict['chr']]= annot_dict[snp].hg19_chr
##                line_list[index_dict['pos']]= annot_dict[snp].hg19_pos
                if not purpose == 'HAPMAP':
                    try:
                        line_list[index_dict['snp']]= annot_dict[snp].lz
                    except KeyError:
                        print("{0} is missing from the dictionary!".format(line_list[index_dict['snp']]))
##                      line_list[index_dict['snp']] = 'ERROR'
                if counter > 5000:
                    #print line_list
                    counter = 0
                counter = 1 + counter
                output.write('\t'.join(line_list)+'\n')
    output.close()
    base, ext = os.path.splitext(table_loc)
    error_loc = base + '_NAMEERRORS.txt'
    efile = open(error_loc,mode="w")
    efile.write('\t'.join(['CHR','RS','cM','POS'])+'\n')
    for error in error_list:
        efile.write('\t'.join(error)+'\n')
    efile.close()
                


def rename_as_necessary(new_orig,ext):
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
    next_rename = new_orig
    if os.path.exists(new_orig + ext):
        print( '''There is already a file named {0}, which means that map_adapt has
already been run on this map file. To avoid overwriting the contents of {0},
the contents of the input file will be written to: '''.format(new_orig + ext))        
        if new_orig[-1].isdigit():
            pieces = new_orig.rpartition('~')
            up1 = int(pieces[2])+ 1
            next_rename = pieces[0]+'~'+ str(up1)
            print(next_rename + ext)
            next_rename = rename_as_necessary(next_rename, ext)
        else:
            next_rename = new_orig + '1'
            print(next_rename + ext)
            next_rename = rename_as_necessary(next_rename, ext)
    return next_rename
    

def read_annot_titles(line):
    title_list = line.strip().split()
    hg18_chr_col = title_list.index(HG18_CHR_TITLE)
    rs_col = title_list.index(RS_ID_TITLE)
    name_col = title_list.index(NAME_TITLE)
    hg18_pos_col = title_list.index(HG18_POS_TITLE)
    band_col = title_list.index(BAND_TITLE)
    hg19_pos_col = title_list.index(HG19_POS_TITLE)
    hg19_chr_col = title_list.index(HG19_CHR_TITLE)
    annot_col = title_list.index(ANNOT_TITLE)
    annot_indices = {'hg18_chr':hg18_chr_col, 'name':name_col,
                    'rs':rs_col,'hg18_pos':hg18_pos_col,
                    'band':band_col,'hg19_pos':hg19_pos_col,
                     'hg19_chr':hg19_chr_col,'annot':annot_col}
    return annot_indices


def locate_annot_dict(build, annot_loc=ANNOT_LOC):
    head, tail = os.path.split(annot_loc)
    dict_loc = os.path.join(head,'annot'+'_'+build+'.dict')
    return dict_loc
    

def build_annot_dict(purpose, annot_dict_loc):
    #global meta
    annot_dict = dict()
    with open(annot_dict_loc, mode="r") as annie:
        for a in annie:
            a_list = a.strip().split()
            ainfo = AnnotInfo._make(a_list)
            
            if ainfo.rs == '0':
                ref = ainfo.name
            if purpose in ['META','FAMILY']:
                annot_dict[ainfo.meta_key]=ainfo
##                annot_dict[ainfo.name]=ainfo
            elif purpose == 'HAPMAP':
                annot_dict[ainfo.name]=ainfo
            elif purpose == 'MAP':
                annot_dict[ainfo.name]=ainfo
            elif purpose == 'LOG':
                annot_dict[ainfo.lz]=ainfo
    return annot_dict

def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global annot_loc, table_loc, meta, build, hapmap, family,fammap
    table_loc = None
    meta = False
    #annot_loc = None
    annot_loc = ANNOT_LOC
    build = BUILD
    hapmap = False
    family = False
    fammap = False
    try: 
        opts, args = getopt.getopt(argv, "h",
                                   ["help","annot=","map=","meta=","build=",
                                    "hapmap","family=","fammap="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("--map"):
            table_loc = arg
        elif opt in ("--meta"):
            table_loc = arg
            meta = True
        elif opt in ("--family"):
            table_loc = arg
            family = True
        elif opt in ("--fammap"):
            table_loc = arg
            fammap = True
        elif opt in ("--hapmap"):
            hapmap = True
        elif opt in ("--annot"):
            annot_loc = None
        elif opt in ("--build"):
            build = arg
        
def usage():
    print('''
USAGE: fix_it.py [FLAG] OBJECT
      FLAG                  DESCRIPTION                         CURRENT DEFAULT
        --meta              filepath of meta table                      
        --map               filepath of plink map file                  
        --annot             filepath of annotation file         {0}                
        --build             genome build ( hg18 or hg19 )       {1}
    -h, --help              display this usage string



'''.format(BUILD, ANNOT_LOC))



def main(argv):
    global annot_loc, table_loc, meta, build, hapmap, purpose, fammap, family
    
    cl_arguments(argv)
    print annot_loc
    if annot_loc is not None:
        print("WRITING ANNOT!")
        write_annot_dict(annot_loc, build)
    if meta:
        purpose = 'META'
    elif hapmap:
        purpose = 'HAPMAP'
    elif family:
        purpose = 'FAMILY'
    elif fammap:
        purpose = 'FAMMAP'
    else:
        purpose = 'MAP'
    dict_loc = locate_annot_dict(build)
    dict_purp = purpose
    if fammap:
        dict_purp = 'FAMILY'
    print purpose
    print dict_purp
    annot_dict = build_annot_dict(dict_purp, dict_loc)
    #annot_dict = build_annot_dict(purpose, dict_loc)
    fix_table(table_loc,annot_dict,purpose,build)

    
if __name__=='__main__':
    main(sys.argv[1:])
