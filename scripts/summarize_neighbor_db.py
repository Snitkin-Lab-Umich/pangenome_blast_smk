import os
import argparse
import pandas as pd

def summarize_neighbor_db(neighbor_db_dir, output_file):
    # for each file in the neighbor_db_dir that ends in .tsv, summarize the number of cases of 2/1/0 neighbors for gene presence/gene absence
    with open(output_file, 'w') as out_f:
        _ = out_f.write('Assembly_Name\tTwo_Neighbors_Present\tOne_Neighbor_Present\tZero_Neighbors_Present\tNot_Found_Present\tTwo_Neighbors_Absent\tOne_Neighbor_Absent\tZero_Neighbors_Absent\tNot_Found_Absent\n')
        for filename in os.listdir(neighbor_db_dir):
            if not filename.endswith('_accessory_neighborfile.tsv'):
                continue
            assembly_name = filename.replace('_accessory_neighborfile.tsv','')
            neighbor_file_path = os.path.join(neighbor_db_dir, filename)
            neighbor_df = pd.read_csv(neighbor_file_path, sep='\t')
            two_n_present, one_n_present, zero_n_present, not_found_present = 0, 0, 0, 0
            two_n_absent, one_n_absent, zero_n_absent, not_found_absent = 0, 0, 0, 0
            for index, row in neighbor_df.iterrows():
                if row['GeneName'] != 'absent':
                    # gene is present
                    neighbors = [row['Neighbor1_GeneName'], row['Neighbor2_GeneName']]
                    n_present = sum([1 for n in neighbors if n != 'no_neighbor'])
                    if neighbors[0] == 'not_found' or neighbors[1] == 'not_found':
                        not_found_present += 1
                    elif n_present == 2:
                        two_n_present += 1
                    elif n_present == 1:
                        one_n_present += 1
                    else:
                        zero_n_present += 1
                else:
                    # gene is absent
                    neighbors = [row['Neighbor1_GeneName'], row['Neighbor2_GeneName']]
                    n_present = sum([1 for n in neighbors if n != 'no_neighbor'])
                    if neighbors[0] == 'not_found' or neighbors[1] == 'not_found':
                        not_found_absent += 1
                    elif n_present == 2:
                        two_n_absent += 1
                    elif n_present == 1:
                        one_n_absent += 1
                    else:
                        zero_n_absent += 1
            out_f.write(f'{assembly_name}\t{two_n_present}\t{one_n_present}\t{zero_n_present}\t{not_found_present}\t{two_n_absent}\t{one_n_absent}\t{zero_n_absent}\t{not_found_absent}\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--neighbor_db','-ndb',type=str,
        help='''Provide a path to a directory containing neighbor data for accessory genes. Files should be in the format [isolate_name]_accessory_neighborfile.tsv.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide a path to an output file.''',
        default=None
        )
    args = parser.parse_args()
    summarize_neighbor_db(args.neighbor_db, args.output)


if __name__ == '__main__':
    main()