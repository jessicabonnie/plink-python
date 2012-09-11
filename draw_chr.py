#! /usr/bin/python2.6
#! ./
'''
Created on May 7, 2012

@author: Jessica Bonnie
'''
import pc_toolbox

TABLE_LOC = '/home/jkb4y/work/data/2012Feb1/eurmeta/eurmeta_06062012_lz_hg19_intersect_noZs.txt'
CHANGE_LOC = '/home/jkb4y/work/data/2012Feb1/eurmeta/chrtable_06062012_hg19_16_intersect.txt'

def edit_table(table_loc, change_loc):
    line1=True
    change_table = open(change_loc, mode="w")
    with open(table_loc, mode = 'r') as table:
        for line in table:
            line_split = line.strip().split()
            if line1:
                index_dict = pc_toolbox.read_meta_titles(line)
                print index_dict
                line1=False
            else:
                if not line_split[index_dict['p']]=='NA':
                    if float(line_split[index_dict['p']]) < 1e-16:
                        if line_split[index_dict['snp']] == 'rs653178':
                            line_split[index_dict['p']] = '1e-16'
                        elif line_split[index_dict['snp']] == 'rs61839660':
                            line_split[index_dict['p']] = '1e-16'
                        elif line_split[index_dict['snp']] == 'rs9273363':
                            line_split[index_dict['p']] = '1e-16'
                        elif line_split[index_dict['snp']] == 'rs2476601':
                            line_split[index_dict['p']] = '1e-16'
                        elif line_split[index_dict['snp']] == 'rs3842727':
                            line_split[index_dict['p']] = '1e-16'
                        elif line_split[index_dict['snp']] == 'rs1701704':
                            line_split[index_dict['p']] = '1.1e-16'
                        ###THESE ARE LOWER
                        elif line_split[index_dict['snp']] == 'rs3087243':
                            line_split[index_dict['p']] = '1e-16'
                        elif line_split[index_dict['snp']] == 'rs12416116':
                            line_split[index_dict['p']] = '1.1e-16'
                        elif line_split[index_dict['snp']] == 'rs3826110':
                            line_split[index_dict['p']] = '1e-16'
                        elif line_split[index_dict['snp']] == 'rs12927355':
                            line_split[index_dict['p']] = '1.1e-16'
                        else:
                            line_split[index_dict['p']] = '1.2e-16'
            change_table.write('\t'.join(line_split) + '\n')
    change_table.close()


def main():
    table_loc = TABLE_LOC
    change_loc = CHANGE_LOC
    edit_table(table_loc, change_loc)

if __name__=='__main__':
    main()
