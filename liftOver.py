#! /ust/bin/python2.6
#! ./
'''
Created May 30, 2012

@author: Jessica Bonnie
'''
import sys
import subprocess

def map2bed(fin, fout):
    fo = open(fout, mode="w")
    fi = open(fin, mode="r")
    for ln in fi:
        chrom, rs, mdist, pos = ln.strip().split()[:4]
        chrom = 'chr' + chrom
        pos = int(pos)
        fo.write('%s\t%d\t%d\t%s\n' % (chrom, pos-1, pos, rs))
    fi.close()
    fo.close()
    return True

def bed2map(fin, fout):
    fo = open(fout, mode="w")
    fi = open(fin,mode="r")
    for ln in fi:
        chrom, pos0, pos1, rs = ln.strip().split()
        chrom = chrom.replace('chr', '')
        fo.write('%s\t%s\t0.0\t%s\n' % (chrom, rs, pos1))
    fo.close()
    fi.close()
    return True

def liftOver(from_build,to_build, bed_input):
    if from_build == 'hg18':
        if to_build == 'hg19':
            chain_file = 'hg18ToHg19.over.chain.gz'
    elif from_build == 'hg19':
        if to_build == 'hg18':
            chain_file = 'hg19ToHg18.over.chain.gz'
        elif to_build == 'hg19':
            chain_file = 'hg19.hg19.all.chain.gz'
    cmd_list = ['liftOver', bed_input, chain_file, '-errorHelp']
    p = subprocess.Popen(cmd_list, bufsize = 0, executable=None,stdin=None,
                         stdout=None,stderr=None, preexec_fn=None,close_fds=False,
                         shell=False,cwd=None,env=None, universal_newlines=False,
                         startupinfo=None, creationflags=0)
    sys.stdout.flush()
    p.wait()
    sys.stdout.flush()

def main(argv):
    from_build = 'hg18'
    to_build = 'hg19'
    map_file = '/home/jkb4y/work/data/HapMap/NoDup/hapmap_removeDup~.bim'
    orig_bed = '/home/jkb4y/work/data/HapMap/NoDup/liftOver/hapmap_removeDup~.bed'
    out_name = '/home/jkb4y/work/data/HapMap/NoDup/liftOver/liftOver'
    out_bed = out_name + '.bed'
    out_map = out_name + '.bim'
    tf = map2bed(map_file, orig_bed)
    liftOver(from_build, to_build, orig_bed)
    #tf = bed2map(out_bed, out_map)




if __name__=='__main__':
    main(sys.argv[1:])
