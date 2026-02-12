import os
import argparse
import pandas as pd

def gene_family_order(input_file, output_file):
    # read in a file of accessory gene families and their neighboring gene families
    # determine the order of the gene families on each scaffold
    # when possible, use the coordinates of the gene family itself
    # when the gene family itself is not present, use the coordinates of the neighboring gene family
    neighbor_df = pd.read_csv(input_file, sep='\t')
    gf_data = {}
    for index, row in neighbor_df.iterrows():
        gf_name = row['GeneFamily']
        gf_scaffold = row['Scaffold']
        if gf_scaffold == 'absent':
            # if both neighbors are present and on the same scaffold, use the average of their coordinates as a rough estimate
            if row['Neighbor1_Scaffold'] != 'absent' and row['Neighbor2_Scaffold'] != 'absent' and row['Neighbor1_Scaffold'] == row['Neighbor2_Scaffold']:
                gf_scaffold = row['Neighbor1_Scaffold']
                gf_genename = f"neighbor_{row['Neighbor1_GeneName']}_{row['Neighbor2_GeneName']}"
                gf_start = (int(row['Neighbor1_Start']) + int(row['Neighbor2_Start'])) // 2
                gf_end = (int(row['Neighbor1_End']) + int(row['Neighbor2_End'])) // 2
            # if only one neighbor is present, or they are on different scaffolds, use the coordinates of the first neighbor
            # (there should never be a case where the first neighbor is absent but the second is present)
            elif row['Neighbor1_Scaffold'] != 'absent':
                gf_scaffold = row['Neighbor1_Scaffold']
                gf_genename = f"neighbor_{row['Neighbor1_GeneName']}"
                gf_start = int(row['Neighbor1_Start']) - 1 
                gf_end = int(row['Neighbor1_End']) - 1
            # if nothing is present, add this to scaffold 'zzz' with a start and end of 0
            else:
                gf_scaffold = 'zzz'
                gf_genename = f"neighbors_absent"
                gf_start = 0
                gf_end = 0
        else:
            gf_genename = row['GeneName']
            gf_start = int(row['Start'])
            gf_end = int(row['End'])
        if gf_scaffold not in gf_data:
            gf_data[gf_scaffold] = []
        gf_data[gf_scaffold].append((gf_name, gf_genename, gf_start, gf_end))
    # sort the gene families on each scaffold by their start position only
    for scaffold in gf_data.keys():
        gf_data[scaffold] = sorted(gf_data[scaffold], key=lambda x: x[2])
    # write the sorted gene families to the output file, with columns for gene family, scaffold, start, and end
    with open(output_file, 'w') as fhout:
        fhout.write('GeneFamily\tScaffold\tStart\tEnd\tGeneName\n')
        for scaffold in sorted(gf_data.keys()):
            for gf_name, gf_genename, gf_start, gf_end in gf_data[scaffold]:
                fhout.write(f'{gf_name}\t{scaffold}\t{gf_start}\t{gf_end}\t{gf_genename}\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a file containing the accessory genes and their neighboring gene families. This should match the output of find_pangenome_neighbors.py.
        Use a chromosome-resolved or long-read assembly when possible.''',
        required=True
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide a path to an output file.''',
        default=None
        )
    args = parser.parse_args()
    gene_family_order(args.input, args.output)


if __name__ == '__main__':
    main()