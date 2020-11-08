import sys, getopt
from align import *


def main(argv):
    referencefile = ''
    queryfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hr:q:o:",["ref_file=","query_file=","output_file="])
    except getopt.GetoptError:
        print('Error! Correct usage:\npython src.py -r <reference_fasta_file> -q <query_fasta_file> [-o <output_csv_file>]')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python src.py -r <reference_fasta_file> -q <query_fasta_file> [-o <output_csv_file>]')
            sys.exit()
        elif opt in ("-r", "--ref_file"):
            referencefile = arg
        elif opt in ("-q", "--query_file"):
            queryfile = arg
        elif opt in ("-o", "--output_file"):
            outputfile = arg
    print('Reference file is ', referencefile)
    print('Query file is ', queryfile)

    all_ref_seqs = parse_fasta(referencefile)
 
    all_seqs = parse_query_file(queryfile)
 
    res_df = align(all_ref_seqs, all_seqs)

    if len(outputfile) == 0:
        output_file = queryfile.split(".")[0] + "_" + referencefile.split("/")[1] + "_results.csv"
    else:
        output_file = outputfile
    print('Output file is ', output_file)
    res_df.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    main(sys.argv[1:])
