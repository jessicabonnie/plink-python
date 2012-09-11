#!/usr/bin/env python
import os
import sys
import re

def in_interval(x,left_endpoint, right_endpoint):
        if left_endpoint <= x and x <= right_endpoint:
                return True
        else:
                return False


gene_region_dict = {}
with open("region_genes_coordinates_NCBI36_wei_min.list","r") as file_in:
        file_in.next()
        for line in file_in:
                split_line = line.strip().split()
                gene_region_dict.update({re.sub("-",".",split_line[4]):[int(split_line[3]), int(split_line[5]), int(split_line[6])]})


nominal_p_cutoff=1e-1
dir_to_look_in = sys.argv[1]
file_out_name = sys.argv[2]
#dir_to_look_in = "Results.B/"
file_list = [filename for filename in os.listdir(dir_to_look_in) if "assoc.linear" in filename and "adjusted" not in filename]
initial_results=[]
for filename in file_list:
        gene_symbol_filename = re.sub("plink\.(.*)\.assoc\.linear","\\1",filename)
        try:
                reg_chrom, reg_start, reg_end = gene_region_dict[gene_symbol_filename]
        except:
                try:
                        reg_chrom, reg_start, reg_end = gene_region_dict[".".join(gene_symbol_filename.split(".")[:-1])]
                except KeyError:
                        continue

        with open(os.path.join(dir_to_look_in,filename), "r") as file_in:
                file_in.next()
                old_P = 1
                old_out=""
                for line in file_in:
                        if "ADD" in line:
                                split_line = line.strip().split()
                                snp_chrom = int(split_line[0])
                                snp_bp = int(split_line[2])
                                if reg_chrom == snp_chrom and in_interval(snp_bp, reg_start, reg_end):
                                        cur_P=float(split_line[8]) 
                                        the_name = split_line[1]
                                        cur_out=[gene_symbol_filename, cur_P, the_name,  snp_chrom, snp_bp, reg_chrom, reg_start, reg_end]
                                        if cur_P < old_P:
                                                old_P = cur_P
                                                old_out = cur_out
                if old_out:
                        initial_results.append(old_out)

                                
print len(set([tuple(item[5:]) for item in initial_results]))
the_dict={}
for item in initial_results:
        key_region = tuple(item[5:])
        if key_region not in the_dict:
                the_dict.update({key_region:item[0:5]})
        else:
                old_p=float(the_dict[key_region][1])
                new_p=float(item[1])
                if new_p < old_p:
                        the_dict.update({key_region:item[0:5]})
with open(file_out_name,"w") as fout:
        fout.write(",".join(["region_chrom","region_start","region_end","gene_regulated","P_value","snp_name","snp_chrom","snp_bp"])+"\n")
        sorted_key =sorted(the_dict.keys())
        for key in sorted_key:
                out=list(key)
                out.extend(the_dict[key])
                out=[str(item) for item in out]
                fout.write(",".join(out)+"\n")

