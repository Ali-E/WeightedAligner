To replicate the reported matching svRNAs of SARS-CoV-2[1] that are orthologs to the ones previously reported for SARS-CoV[2], run the command:

python src.py -r Data/coronavirus_MT327745.fa -q Data/three_svRNA.fa -o match_results.csv

The default scores used in alignment are as follows:

Match: +1, indel: -1, mismatch (other than AG and CT): -1, A/G mismatch: -0.5, C/T mismatch: -0.75

To change these scores use -a (for A/G mismatch) and -c (for C/T mismatch), for example:


python src.py -r Data/coronavirus_MT327745.fa -q Data/three_svRNA.fa -o match_results.csv -a 0.0 -c -1



If you use the code of this aligner (or parts of it) please cite [1].

[1] Boroojeny, A. E., & Chitsaz, H. (2020). SARS-CoV-2 orthologs of pathogenesis-involved small viral RNAs of SARS-CoV. arXiv preprint arXiv:2007.05859.

[2] Morales, L., Oliveros, J. C., Fernandez-Delgado, R., Robert tenOever, B., Enjuanes, L., & Sola, I. (2017). SARS-CoV-encoded small RNAs contribute to infection-associated lung pathology. Cell host & microbe, 21(3), 344-355.
