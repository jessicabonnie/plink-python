#! /ust/bin/python2.6
#! ./
'''
Created Jan 25, 2012

@author: Jessica Bonnie
'''
import os
import sys
import pc_toolbox
import fix_it
import achilleas_yank

P_TITLE = 'P'


def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global tablefold, superfolder, cell_type
    tablefold = None
    superfolder = None
    cell_type = None
    try: 
        opts, args = getopt.getopt(argv, "h",
                                   ["help","tablefolder=","superfolder=",
                                    "cell-type=","build=",
                                    "hapmap","family=","fammap="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("--outfolder"):
            superfolder = arg
        elif opt in ("--tablefolder"):
            tablefold = arg
        elif opt in ("--cell-type"):
            cell_type = arg




def write_prelude(prelude_loc, gene):
    base, ext = os.path.splitext(prelude_loc)
    prelude_gene_loc = base + '_'+gene + ".R"
    prelude = open(prelude_gene_loc, mode = "w")
    prelude.write("gene = "+gene+"\n")
    gen = open(prelude_loc, mode = "r")
    for line in gen:
        prelude.write(line)
    return prelude_gene_loc
    






def main(argv):
    global tablefold, superfolder, cell_type
    #cl_arguments(argv)
    cell_type = 'B'
    region_list = '/home/jkb4y/work/data/Region_Lists/hg19/T1D_regions_chr17_18_19.txt'#hg19_dbase_and_sig.txt'
    superfolder = '/home/jkb4y/ubs/work/results/Achilleas/hg19/'
    tablefold = '/home/jkb4y/ubs/work/data/Achilleas/cis-eQTL/'
    table_loc = os.path.join(tablefold,'{0}_cis_eqtls_permwBP_JB.txt'.format(cell_type))
    cellfolder = os.path.join(superfolder, cell_type)
    outfolder = os.path.join(cellfolder, 'RegionGraphs')

    draw_clist= ['python', 'illustrate.py', '--build', 'hg19',
                 '--region-list', region_list,
                 '--eqtl',cell_type,'--outfolder', outfolder,
                 '--table', table_loc, '--max-bp','0.05']
    print draw_clist
    cmd = ' '.join(draw_clist)
    os.system(cmd)




if __name__=='__main__':
    main(sys.argv[1:])
