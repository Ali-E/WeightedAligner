#!/bin/bash
# rm coronavirus_final_results_filenames_new.txt
# rm coronavirus_filenames_final.txt

# for filename in coronavirus_sequences_new/*_results.csv
# for filename in coronavirus_sequences_new/*.fa
for filename in RNA-seq/*.fa
do
    # echo "$filename" >> coronavirus_final_results_filenames_new.txt

    # echo "$filename" >> coronavirus_filenames_final.txt
    echo "$filename" >> coronavirus_RNA_seq.txt
    python align.py $filename three_svRNA_orthologs_final_truncated.fa 
done

