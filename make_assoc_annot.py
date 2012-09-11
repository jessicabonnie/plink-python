import os
import pc_toolbox

def main():
    build = 'hg19'
    table_loc = '/home/jkb4y/ubs/work/results/July2012/Logistic_pc8/AA_8.assoc.logistic'
    
    annot_dict = pc_toolbox.create_annot_dict(build,'LOG')
    base, ext = os.path.splitext(table_loc)
    base2, ext2 = os.path.splitext(base)
    new_loc = base2 + '_lz'+ext2+ext
    index_dict, line = pc_toolbox.create_index_dict(table_loc,'assoc')
    pc_toolbox.give_table_annotation(table_loc,new_loc,index_dict,annot_dict)


if __name__ == '__main__':
    main()
