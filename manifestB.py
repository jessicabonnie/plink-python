#! /ust/bin/python2.7
#! ./
'''
Created May 25, 2012

@author: Jessica Bonnie
'''
from collections import namedtuple
import fix_it
import shutil
import os
import sys
import getopt
MANIFEST_LOC = '/home/jkb4y/work/data/Immuno_BeadChip_11419691_Manifest_B.txt'

CHR_TITLE = 'Chr'
IM_TITLE = 'Name'
POS_TITLE = 'MapInfo'
BUILD_TITLE = 'GenomeBuild'

BUILD = 'hg18'

AnnotInfo = fix_it.AnnotInfo
ManifestInfo = namedtuple("ManifestInfo","name,build,chro,pos")
TABLE_LOC = '/home/jkb4y/work/data/HapMap/NoDupSAMP/hapmap_removeDup.bim'

def read_manifest(manifest_indices, line):
    mlist = line.strip().split('\t')
    chro=mlist[manifest_indices['chr']].replace('X','23').replace('Y','24').replace('XY','23')
    name=mlist[manifest_indices['name']]
    pos=mlist[manifest_indices['pos']]
    gbuild=mlist[manifest_indices['build']]
    if 36 <= float(gbuild) < 37:
        hgbuild = 'hg18'
    elif 37 <= float(gbuild) < 38:
    #elif gbuild == '37':
        hgbuild = 'hg19'
    else:
        hgbuild = 'ERROR'
    manifest_info = ManifestInfo(name=name,build=hgbuild,chro=chro,pos=pos)
    return manifest_info



def write_manifest_dict(manifest_loc):
    dict_loc = locate_manifest_dict(manifest_loc)
    manifest_dict = open(dict_loc, mode="w")
    with open(manifest_loc, mode = "r") as manifest:
        line1 = True
        for mline in manifest:
            if line1:
                m_indices = read_manifest_titles(mline)
                line1 = False
            else:
                m_info = read_manifest(m_indices,mline)
                manifest_dict.write('\t'.join(m_info) +'\n')
    manifest_dict.close()

def read_manifest_titles(line):
    title_list = line.strip().split('\t')
    chr_col = title_list.index(CHR_TITLE)
    name_col = title_list.index(IM_TITLE)
    pos_col = title_list.index(POS_TITLE)
    build_col = title_list.index(BUILD_TITLE)
    manifest_indices = {'chr':chr_col, 'name':name_col,
                        'pos':pos_col,'build':build_col}
    return manifest_indices

def locate_manifest_dict(manifest_loc=MANIFEST_LOC):
    head, tail = os.path.split(manifest_loc)
    dict_loc = os.path.join(head,'manifest_B'+'.dict')
    return dict_loc

def locate_table_key(table_loc):
    base, ext = os.path.splitext(table_loc)
    key_loc = base + '_key.tbl'
    return key_loc
    

def build_manifest_dict(manifest_dict_loc):
    #global meta
    manifest_dict = dict()
    with open(manifest_dict_loc, mode="r") as mannie:
        for m in mannie:
            m_list = m.strip().split()
            minfo = ManifestInfo._make(m_list)
            ref = minfo.name
            manifest_dict[ref]=minfo
    return manifest_dict



def fix_table(table_loc, build, manifest_dict, annot_dict):

    (basepath, ext) = os.path.splitext(table_loc)
    orig_rename = basepath + '~'
    new_loc = str(rename_as_necessary(orig_rename, ext)) + ext
    output_loc = new_loc


    with open(table_loc, mode="r") as table:
        output = open(output_loc, mode="w")
        index_dict = {'chr':0,'pos':3,'snp':1}
        error_list = list()
        for line in table:
            line_list = line.strip().split()
            snp = line_list[index_dict['snp']]
            if manifest_dict[snp].build == build:
                line_list[index_dict['chr']]= manifest_dict[snp].chro.replace('X','23').replace('Y','24').replace('XY','23')
                line_list[index_dict['pos']]= manifest_dict[snp].pos
            elif snp in ['rs9279834','rs10636228','rs17510415',
                         'rs17539100','rs17510360','rs5743293','rs35768562']:
                line_list[index_dict['pos']]= str(int(line_list[index_dict['pos']])+1)
            else:
                try:
                    if build == 'hg18':
                        line_list[index_dict['chr']]= annot_dict[snp].hg18_chr.replace('X','23').replace('Y','24').replace('XY','23')
                        line_list[index_dict['pos']]= annot_dict[snp].hg18_pos
                    else:
                        line_list[index_dict['chr']]= annot_dict[snp].hg19_chr.replace('X','23').replace('Y','24').replace('XY','23')
                        line_list[index_dict['pos']]= annot_dict[snp].hg19_pos
                except KeyError:
                    error_list.append(list(line_list))
            output.write('\t'.join(line_list)+'\n')
    output.close()
    base, ext = os.path.splitext(new_loc)
    error_loc = base + '_NAMEERRORS.txt'
    efile = open(error_loc,mode="w")
    efile.write('\t'.join(['CHR','RS','cM','POS'])+'\n')
    for error in error_list:
        efile.write('\t'.join(error)+'\n')
    efile.close()
                

def create_key(table_loc, manifest_dict, annot_dict, key_loc):
    KeyInfo = namedtuple("KeyInfo","im,rs,hg18_chr,hg18_pos,hg19_chr,hg19_pos,mbuild")
    key = open(key_loc, mode="w")
    title_list = ['imchip_id','rs_id','hg18_chr','hg18_pos','hg19_chr','hg19_pos','manifest_build']
    key.write('\t'.join(title_list)+'\n')
    with open(table_loc, mode="r") as table:
        index_dict = {'chr':0,'pos':3,'snp':1}
        error_list = list()
        
        for line in table:
            
            line_list = line.strip().split()
            im = line_list[index_dict['snp']]
            mbuild = manifest_dict[im].build
            rs = 'NA'
            hg18_chr = 'NA'
            hg19_chr = 'NA'
            hg18_pos = 'NA'
            hg19_pos = 'NA'
            print mbuild
            if mbuild == 'hg18':
                hg18_chr = manifest_dict[im].chro
                hg18_pos = manifest_dict[im].pos
                try:
                    hg19_chr = annot_dict[im].hg19_chr
                    hg19_pos = annot_dict[im].hg19_pos
                    rs = annot_dict[im].rs
                    
                except KeyError:
                    rs = 'NA'
                    hg19_chr = 'NA'
                    hg19_pos = 'NA'
##                    error_list.append(list(line_list))
##                    error_list.extend(mbuild)
                    line_list.append(mbuild)
                    error_list.append(list(line_list))
            elif im in ['rs9279834','rs10636228','rs17510415',
                        'rs17539100','rs17510360','rs5743293','rs35768562']:
                hg18_chr = line_list[index_dict['chr']]
                hg18_pos = str(int(line_list[index_dict['pos']])+1)
                hg19_chr = 'NA'
                hg19_pos = 'NA'
                rs = im
            elif mbuild == 'hg19':
                hg19_chr = manifest_dict[im].chro
                hg19_pos = manifest_dict[im].pos
                try:
                    hg18_chr = annot_dict[im].hg18_chr
                    hg18_pos = annot_dict[im].hg18_pos
                    rs = annot_dict[im].rs
                except KeyError:
                    line_list.append(mbuild)
                    error_list.append(list(line_list))
                    rs = 'NA'
                    hg18_chr = 'NA'
                    hg18_pos = 'NA'
            keyinfo = KeyInfo(im=im,rs=rs,hg18_chr=hg18_chr,
                              hg18_pos=hg18_pos, hg19_chr=hg19_chr,
                              hg19_pos=hg19_pos,mbuild=mbuild)

            key.write('\t'.join(keyinfo)+'\n')
    key.close()
    base, ext = os.path.splitext(table_loc)
    error_loc = base + '_KEYERRORS.txt'
    efile = open(error_loc,mode="w")
    efile.write('\t'.join(['CHR','RS','cM','POS','A1','A2','manifest_build'])+'\n')
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
the contents of the new output file will be written to: '''.format(new_orig + ext))        
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
    






def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global manifest_loc, build, table_loc
    manifest_loc = None
    #manifest_loc = MANIFEST_LOC
    build = BUILD
    table_loc = TABLE_LOC
    try: 
        opts, args = getopt.getopt(argv, "h",
                                   ["help","manifest=","build=","table="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("--table"):
            table_loc = arg
        elif opt in ("--manifest"):
            manifest_loc = arg
        elif opt in ("--build"):
            build = arg





def main(argv):
    global manifest_loc, table_loc, build
    cl_arguments(argv)
    if manifest_loc is not None:
        print("WRITING MANIFEST DICT!")
        write_manifest_dict(manifest_loc)
    dict_loc = locate_manifest_dict()
    
    manifest_dict = build_manifest_dict(dict_loc)
    annot_dict_loc = fix_it.locate_annot_dict(build)
    annot_dict = fix_it.build_annot_dict('HAPMAP', annot_dict_loc)
    fix_table(table_loc, build, manifest_dict, annot_dict)
    key_loc = locate_table_key(table_loc)
    create_key(table_loc, manifest_dict, annot_dict, key_loc)
    
if __name__=='__main__':
    main(sys.argv[1:])
