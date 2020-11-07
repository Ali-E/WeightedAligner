import pandas as pd
import sys
import numpy as np



if __name__ == "__main__":
    filename_file = sys.argv[1]
    
    filenames = []
    with open(filename_file) as f:
        lines = f.readlines()
        filenames = [line[:-1] for line in lines]

    svRNA1_data = []
    svRNA2_data = []
    svRNA3_data = []
    svRNA_data = []
    bad_files = []
    missing_files = []

    for filename in filenames:
        new_df = pd.read_csv(filename, delimiter='\t')
        print(new_df)
       
        score_col = 1
        if len(new_df):
            svRNA1_data.append(new_df.iloc[0,score_col])
            svRNA2_data.append(new_df.iloc[1,score_col])
            svRNA3_data.append(new_df.iloc[2,score_col])
            if new_df.iloc[1,score_col] == 1:
                svRNA_data.append((filename.split('/')[1].split('_')[1].split('.')[0], new_df.iloc[0,6] - new_df.iloc[0,2] + 1, new_df.iloc[1,6] - new_df.iloc[1,2] + 1, new_df.iloc[2,6] - new_df.iloc[2,2] + 1)) # without ref_name column
                # svRNA_data.append((filename.split('/')[1].split('_')[5].split('.')[0], new_df.iloc[0,7] - new_df.iloc[0,3] + 1, new_df.iloc[1,7] - new_df.iloc[1,3] + 1, new_df.iloc[2,7] - new_df.iloc[2,3] + 1)) # with ref_name column


            print(new_df.iloc[0,score_col])
            print(new_df.iloc[1,score_col])
            print(new_df.iloc[2,score_col])
            if new_df.iloc[1,score_col] != 1:
                # print(filename)
                bad_files.append(filename)
                # exit(0)
        else:
            missing_files.append(filename)


    print(np.mean(svRNA1_data))
    print(np.mean(svRNA2_data))
    print(np.mean(svRNA3_data))
    print(len(svRNA1_data))
            

    print('bad: ', bad_files)
    print('missing: ', missing_files)
    for missing in missing_files:
        print(missing)


    df_new = pd.DataFrame(svRNA_data, columns=['NCBI ID', 'svRNA1', 'svRNA2', 'svRNA3'])
    # df_new.to_csv('table_svRNA_orthologs_final_all_refs_new.csv', sep='\t')


