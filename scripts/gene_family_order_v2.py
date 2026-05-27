import os
import argparse
import pandas as pd
import gffutils as gff

def gene_family_order(input_file, gff_file, output_file, splice_mode=False):
    # read in a file of accessory gene families and their neighboring gene families
    # determine the order of the gene families on each scaffold
    # when possible, use the coordinates of the gene family itself
    # when the gene family itself is not present, use the coordinates of the neighboring gene family
    # all coordinates should be extracted from the provided gff file, rather than the coordinates present in the input neighbor file
    # this gives me the freedom to use either spliced or unspliced coordinates as needed
    neighbor_df = pd.read_csv(input_file, sep='\t')
    gff_db = gff.create_db(gff_file,dbfn=":memory:",force=True,keep_order=False,merge_strategy="create_unique",sort_attribute_values=True,from_string=False)
    gf_data = {}
    for index, row in neighbor_df.iterrows():
        gf_name = row['GeneFamily']
        gf_scaffold = row['Scaffold']
        if gf_scaffold == 'absent':
            # if both neighbors are present and on the same scaffold, use the average of their coordinates as a rough estimate
            if row['Neighbor1_Scaffold'] != 'absent' and row['Neighbor2_Scaffold'] != 'absent' and row['Neighbor1_Scaffold'] == row['Neighbor2_Scaffold']:
                gf_scaffold = row['Neighbor1_Scaffold']
                gf_genename = f"neighbor_{row['Neighbor1_GeneName']}_{row['Neighbor2_GeneName']}"
                neighbor1_scaffold, neighbor1_start, neighbor1_end = get_gene_data_from_gff(gff_db, row['Neighbor1_GeneName'], splice_mode=splice_mode)
                neighbor2_scaffold, neighbor2_start, neighbor2_end = get_gene_data_from_gff(gff_db, row['Neighbor2_GeneName'], splice_mode=splice_mode)
                # check to make sure the correct annotations were pulled, based on the scaffold names
                if neighbor1_scaffold != gf_scaffold or neighbor2_scaffold != gf_scaffold:
                    print(f'Error: gff returned {neighbor1_scaffold} and {neighbor2_scaffold} when {gf_scaffold} was expected.')
                    quit(1)
                gf_start = min([int(neighbor1_start),int(neighbor2_start)])
                gf_end = max([int(neighbor1_end),int(neighbor2_end)])
            # if only one neighbor is present, or they are on different scaffolds, use the coordinates of the first neighbor
            # (there should never be a case where the first neighbor is absent but the second is present)
            elif row['Neighbor1_Scaffold'] != 'absent':
                gf_scaffold = row['Neighbor1_Scaffold']
                gf_genename = f"neighbor_{row['Neighbor1_GeneName']}"
                neighbor1_scaffold, neighbor1_start, neighbor1_end = get_gene_data_from_gff(gff_db, row['Neighbor1_GeneName'], splice_mode=splice_mode)
                if neighbor1_scaffold != gf_scaffold:
                    print(f'Error: gff returned {neighbor1_scaffold} when {gf_scaffold} was expected.')
                    quit(1)
                gf_start = int(neighbor1_start) 
                gf_end = int(neighbor1_end)
            # if nothing is present, add this to scaffold 'zzz' with a start and end of 0
            else:
                gf_scaffold = 'zzz'
                gf_genename = f"neighbors_absent"
                gf_start = 0
                gf_end = 0
        # if the gene family itself is present, use its coordinates
        else:
            gf_genename = row['GeneName']
            gene_scaffold, gene_start, gene_end = get_gene_data_from_gff(gff_db, gf_genename, splice_mode=splice_mode)
            if gene_scaffold != gf_scaffold:
                print(f'Error: gff returned {gene_scaffold} when {gf_scaffold} was expected.')
                quit(1)
            gf_start = int(gene_start)
            gf_end = int(gene_end)
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

def get_gene_data_from_gff(gff_db, gene_name_cds, splice_mode=False):
    # take a gffutils database and a gene name corresponding to a cds feature
    # return the scaffold, start, and end coordinates of the mRNA feature corresponding to this gene name
    # if splice_mode is True, instead return the coordinates of the CDS feature itself (since this assumes a spliced gff was used)
    gene_feature_cds = gff_db[gene_name_cds]
    if splice_mode:
        return gene_feature_cds.seqid, gene_feature_cds.start, gene_feature_cds.end
    else:
        gene_feature_mrna = gff_db[gene_feature_cds.attributes['Parent'][0]]
        return gene_feature_mrna.seqid, gene_feature_mrna.start, gene_feature_mrna.end

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
    parser.add_argument(
        '--gff','-g',type=str,
        help='''Provide a path to a GFF file for the assembly used in the input file. This is used to get the coordinates of gene families that are not in the input file.''',
        required=True
        )
    parser.add_argument(
        '--splice_mode','-s',action='store_true',
        help='''If set, use spliced coordinates instead of unspliced coordinates for gene families that are in the input file.''',
        default=False
        )
    args = parser.parse_args()
    gene_family_order(args.input, args.gff, args.output, splice_mode=args.splice_mode)


if __name__ == '__main__':
    main()