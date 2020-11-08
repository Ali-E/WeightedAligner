import sys, getopt
from align import *


def main(argv):
    referencefile = ''
    queryfile = ''
    outputfile = ''
    AG_score = -0.5
    CT_score = -0.75
    try:
        opts, args = getopt.getopt(argv,"hr:q:o:a:c:",["ref_file=","query_file=","output_file=","AG_score=","CT_score="])
    except getopt.GetoptError:
        print('Error! Correct usage:\npython src.py -r <reference_fasta_file> -q <query_fasta_file> [-o <output_csv_file> -a <score for AG mismatch> -c <score for CT mismatch>]')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python src.py -r <reference_fasta_file> -q <query_fasta_file> [-o <output_csv_file> -a <score for AG mismatch> -c <score for CT mismatch>]')
            sys.exit()
        elif opt in ("-r", "--ref_file"):
            referencefile = arg
        elif opt in ("-q", "--query_file"):
            queryfile = arg
        elif opt in ("-o", "--output_file"):
            outputfile = arg
        elif opt in ("-a", "--AG_score"):
            AG_score = float(arg)
        elif opt in ("-c", "--CT_score"):
            CT_score = float(arg)
    print('Reference file is ', referencefile)
    print('Query file is ', queryfile)

    print('AG score: ', AG_score)
    print('CT score: ', CT_score)

    all_ref_seqs = parse_fasta(referencefile)
 
    all_seqs = parse_query_file(queryfile)
 
    res_df = align(all_ref_seqs, all_seqs, AG=AG_score, CT=CT_score)

    if len(outputfile) == 0:
        output_file = queryfile.split(".")[0] + "_" + referencefile.split("/")[1] + "_results.csv"
    else:
        output_file = outputfile
    print('Output file is ', output_file)
    res_df.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    main(sys.argv[1:])
