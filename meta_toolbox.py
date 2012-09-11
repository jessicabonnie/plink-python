#! /ust/bin/python2.7
#! ./
'''
Created Jan 25, 2012

@author: Jessica Bonnie
'''
from collections import namedtuple
import pc_toolbox
from operator import itemgetter


AnnotInfo = namedtuple("AnnotInfo","chro,pos,hg18,rs,lz,band")
ANNOT_LOC = '/home/jkb4y/work/data/quinlan-immunochip-snps-annotated-2011-Dec-15_edit.txt'
ANNOT_DICT = '/home/jkb4y/work/data/annotation_lz.dict'

META_CHR_TITLE = 'CHR'
SNP_TITLE = 'MarkerName'
P_TITLE = 'P-value'
POSITION_TITLE = 'POS'
Z_TITLE = 'Zscore'
WEIGHT_TITLE = 'Weight'
A1_TITLE = 'Allele1'
A2_TITLE = 'Allele2'

HG18_CHR_TITLE = 'hg18_chrom'
HG18_ID_TITLE = 'name'
RS_ID_TITLE = 'rsID'
HG18_POS_TITLE = 'hg18_end'
BAND_TITLE = 'band'

def read_meta_titles(line):
    title_list = line.strip().split()
    meta_chr_col = title_list.index(META_CHR_TITLE)
    p_col = title_list.index(P_TITLE)
    snp_col = title_list.index(SNP_TITLE)
    position_col = title_list.index(POSITION_TITLE)
    z_col = title_list.index(Z_TITLE)
    weight = title_list.index(WEIGHT_TITLE)
    meta_a = title_list.index(A1_TITLE)
    a2 = title_list.index(A2_TITLE)
    meta_indices = {'chr':meta_chr_col, 'p':p_col,
                    'snp':snp_col,'z':z_col,
                    'pos':position_col,'weight':weight,
                    'a1':meta_a, 'a2':a2,'pvalcol':P_TITLE,'markercol':SNP_TITLE}
##    return meta_indices
##def read_meta_titles(line):
##    title_list = line.strip().split()
##    meta_chr_col = title_list.index(META_CHR_TITLE)
##    p_col = title_list.index(P_TITLE)
##    snp_col = title_list.index(SNP_TITLE)
##    position_col = title_list.index(POSITION_TITLE)
##    z_col = title_list.index(Z_TITLE)
##    weight = title_list.index(WEIGHT_TITLE)
##    meta_a = title_list.index(A1_TITLE)
##    a2 = title_list.index(A2_TITLE)
##    meta_indices = {'chr':meta_chr_col, 'p':p_col,
##                    'snp':snp_col,'z':z_col,
##                    'pos':position_col,'weight':weight,
##                    'a1':meta_a, 'a2':a2}
    print("Index Dict:")
    print meta_indices
    return meta_indices


def meta_tuple(meta_line, index_dict, annot_dict=None):
    '''creates named tuple of the form (chro=chromosome,
                                        p=p-value,
                                        pos=position
                                        snp=SNP,
                                        lz=Locuszoom SNP
                                        z=z-score,
                                        a1=A1 allele,
                                        a2=A2 allele,
                                        w= weight,
                                        abs_z=|z-score|)

    Keyword arguments:
    meta_line -- line from meta_table

    Returns - namedtuple

    '''
    line_list = meta_line.strip().split()
    p_store = line_list[index_dict['p']]
    if pc_toolbox.is_number(p_store):
        p_store = float(p_store)
    im_snp, lz_snp = correct_snp(index_dict,line_list, annot_dict)
    if annot_dict is None:
        band = '--'
    else:
        band = annot_dict[im_snp].band
    #annot_i = getIndexOfTuple(annot_list,0,im_snp)
    abs_z = abs(float(line_list[index_dict['z']]))
    LineInfo = namedtuple('LineInfo', 'chro,p,pos,snp,lz,z,a1,a2,w,abs_z,band')
    meta_tuple = LineInfo(chro=line_list[index_dict['chr']],
                          p=line_list[index_dict['p']],
##                          p=p_store,
                          pos=line_list[index_dict['pos']],
                          snp=im_snp,
                          lz=lz_snp,
                          z=line_list[index_dict['z']],
                          abs_z=abs_z,
                          a1=line_list[index_dict['a1']],
                          a2=line_list[index_dict['a2']],
                          w=line_list[index_dict['weight']],
                          band=band)
    return meta_tuple

def read_annot_titles(line):
    title_list = line.strip().split()
    hg18_chr_col = title_list.index(HG18_CHR_TITLE)
    rs_col = title_list.index(RS_ID_TITLE)
    hg18_col = title_list.index(HG18_ID_TITLE)
    pos_col = title_list.index(HG18_POS_TITLE)
    band_col = title_list.index(BAND_TITLE)
    annot_indices = {'chro':hg18_chr_col, 'hg18':hg18_col,
                    'rs':rs_col,'pos':pos_col,
                    'band':band_col}
    return annot_indices

def read_annot(annot_loc=ANNOT_LOC):
    counter = 0
##    new_annot = '/home/jkb4y/work/data/imchip_annotation_short.txt'
##    new = open(new_annot, mode = "w")
##    titles = [HG18_CHR_TITLE,HG18_POS_TITLE,HG18_ID_TITLE,RS_ID_TITLE","lz_id","band"]
##    new.write('\t'.join(titles)+"\n")
    with open(annot_loc, mode = "r") as annot:
        line1 = True
        annot_dict=dict()
        annot_list = list()
        keep_rs_list = list()
        AnnotInfo = namedtuple("AnnotInfo","chro,pos,hg18,rs,lz,band")
        for aline in annot:
            if line1:
                a_indices = read_annot_titles(aline)
                line1 = False
            else:
                alist = aline.strip().split()
                hg18=alist[a_indices['hg18']]
                rs=alist[a_indices['rs']]
                pos=alist[a_indices['pos']]
                band=alist[a_indices['band']]
                chro=alist[a_indices['chro']]
                ref = rs
                if rs == '0':
                    ref = hg18
##                if hg18 == rs:
##                    lz=hg18
##                    keep_rs_list.append(hg18)
                if hg18.startswith("rs"):
                    lz=hg18
                    keep_rs_list.append(hg18)
                else:
                    lz = "chr"+chro+":"+pos
                annot_info = AnnotInfo(chro=chro,pos=pos,
                                       hg18=hg18,
                                       rs=rs,
                                       lz=lz,band=band)
##                if not ref in keep_rs_list:
                annot_dict[ref]=annot_info
                annot_dict[hg18]=annot_info
                if counter > 7000:
                    print annot_dict[ref]
                    counter = 0
                counter = counter + 1
##                new.write("\t".join([chro,pos,hg18,rs,lz,band])+"\n")
##    new.close()
    return annot_dict

def getIndexOfTuple(listy, index, value):
    for pos,t in enumerate(listy):
        if t[index] == value:
            return pos
    # Matches behavior of list.index
    raise ValueError("list.index(x): x not in list at value={0}".format(value))

def read_annot_dict(annot_dict_loc=ANNOT_DICT):
    annot_dict = dict()
    with open(annot_dict_loc, mode="r") as annie:
        for a in annie:
            a_list = a.strip().split()
            ainfo = AnnotInfo._make(a_list)
            ref = ainfo.rs
            if ainfo.rs == '0':
                ref = ainfo.hg18
            annot_dict[ref]=ainfo
            annot_dict[ainfo.hg18]=ainfo
##            if not ref == ainfo.hg18:
##                annot_dict[ainfo.hg18]=ainfo
    return annot_dict

def correct_snp(index_dict, line_list, annot_dict):
    snp = line_list[index_dict['snp']]
    if annot_dict is None:
        if snp.startswith("rs"):
            return snp, snp
        else:
            lz_snp = 'chr'+line_list[index_dict['chr']]+':'+line_list[index_dict['pos']]
            return snp, lz_snp
    else:
        return snp, annot_dict[snp].lz
##    if snp in keep_rs_list:
##        return snp, snp
##    return snp, annot_dict[snp].lz
##    if snp in keep_rs_list:
##        if counter > 3:
##            print("{0} is in keep list!".format(snp))
##            counter = 0
##        counter = counter + 1
##        return snp , snp
##    else:
##        if counter > 3:
##            print("{0} is not in keep list!".format(snp))
##            counter = 0
##        counter = counter + 1
##        annot_i = getIndexOfTuple(annot_list,1,snp)
##        lz_snp = 'chr'+annot_list[annot_i].chro+':'+annot_list[annot_i].pos
##        return snp, lz_snp
##    pos = line_list[index_dict['pos']]
##    ch = line_list[index_dict['chr']]
    #if snp in repair_dict
##    if snp.startswith('rs'):
##        return snp , snp
##    else:
##        lz_snp = 'chr'+ch+':'+pos
##        return snp, lz_snp
