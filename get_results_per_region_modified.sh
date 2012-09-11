#!/bin/bash
for item in {CD4,CD8,B,MONO,NK};do
        echo $item
        python my_filter_per_region_modified.py Results.${item} ${item}_results_per_region_modified.csv
done
