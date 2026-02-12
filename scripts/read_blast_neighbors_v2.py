
import argparse
import pandas as pd
from itertools import filterfalse

def write_output_v2(outdata, output_file):
    with open(output_file, 'w') as out_f:
        _ = out_f.write('Query_SeqID\tSubject_SeqID\tExpected_Neighbors\tMax_Neighbors_in_BLAST\n')
        for query_seqid, subject_hit_data_list in outdata.items():
            if subject_hit_data_list == []:
                subject_hit_data_list = [('None','0','0')]
            for subject_hit_data in subject_hit_data_list:
                assembly_name, neighbors_expected, max_neighbors_of_blast_hit = subject_hit_data
                out_f.write(f'{query_seqid}\t{assembly_name}\t{neighbors_expected}\t{max_neighbors_of_blast_hit}\n')

def evaluate_blast_hits_v2(blastdf, neighbor_db_dir, minimum_identity=0.8, minimum_evalue=1e-5, minimum_coverage=0.8, 
    max_distance=50000):
    # for each query sequence, return the list of subject sequence IDs that match the provided criteria
    # take all hits for a specific subject sequence ID into account, as long as they are above the thresholds
    outdata = {}
    for query_seqid, query_blastdf in blastdf.groupby('query_seqid'):
        #print(f'Evaluating hits for query {query_seqid}')
        subject_hit_data_list = find_neighboring_best_hit_v2(query_blastdf, minimum_identity, minimum_evalue, minimum_coverage, max_distance, neighbor_db_dir)
        outdata[query_seqid] = subject_hit_data_list
    return outdata

def find_neighboring_best_hit_v2(query_blastdf, minimum_identity, minimum_evalue, minimum_coverage, max_distance, neighbor_db_dir):
    # use a dictionary of neighboring gene families to identify if there is a high-quality hit in the expected location
    subject_hit_data = []
    # extract all hits that meet the criteria for each assembly
    query_blastdf['assembly_name'] = query_blastdf['subject_seqid'].apply(lambda x: x.split('.scaffolds')[0])
    query_blastdf['subject_scaffold'] = query_blastdf['subject_seqid'].apply(lambda x: x.split('.scaffolds_')[1])
    for assembly_name, assembly_blastdf in query_blastdf.groupby('assembly_name'):
        assembly_blastdf['coverage'] = assembly_blastdf['alignment_length'] / assembly_blastdf['query_length']
        blastdf_filtered = assembly_blastdf[
            (assembly_blastdf['percent_identity'] >= minimum_identity * 100) &
            (assembly_blastdf['evalue'] <= minimum_evalue) &
            (assembly_blastdf['coverage'] >= minimum_coverage)
        ]
        if blastdf_filtered.empty:
            continue
        neighbor_df = pd.read_csv(f'{neighbor_db_dir}/{assembly_name}_accessory_neighborfile.tsv', sep='\t')
        neighbors_expected, max_neighbors_of_blast_hit = gff_neighbor_search_v2(blastdf_filtered, neighbor_df, max_distance)
        subject_hit_data.append((assembly_name, neighbors_expected, max_neighbors_of_blast_hit))
    return subject_hit_data


def gff_neighbor_search_v2(blastdf, neighbor_df, max_distance):
    # blastdf should consist of a single gene family query and a single assembly name
    # first, extract the gene family and the assembly name
    query_gene_fam = blastdf['query_seqid'].iloc[0]
    subject_assembly = blastdf['assembly_name'].iloc[0]
    neighbors_expected = 0
    #print(f'Checking neighbors for gene family {query_gene_fam} in assembly {subject_assembly}')
    # get the corresponding row from the neighbor_df
    neighbor_data = neighbor_df[neighbor_df['GeneFamily'] == query_gene_fam]
    if neighbor_data['Neighbor1_Scaffold'].iloc[0] == 'absent' and neighbor_data['Neighbor2_Scaffold'].iloc[0] != 'absent':
        print('unexpected neighbor data format, neighbor 1 is absent but neighbor 2 is present!')
        quit(1)
    # count number of non-absent neighbors
    if neighbor_data['Neighbor1_Scaffold'].iloc[0] == 'absent':
        neighbors_expected = 0
    else:
        neighbors_expected = 1
        if neighbor_data['Neighbor2_Scaffold'].iloc[0] != 'absent':
            neighbors_expected = 2
    if neighbors_expected == 0:
        # do not search further if there are no expected neighbors
        return(neighbors_expected, 0)
    # for each neighbor position, check if it is on the same scaffold and within the max_distance of any of the BLAST hits
    # return the number of neighbors that meet this requirement
    blast_results = []
    for _, row in blastdf.iterrows():
        blast_outcome = 0
        if row['subject_scaffold'] == neighbor_data['Neighbor1_Scaffold'].iloc[0]:
            # check distance
            distance = min(abs(row['subject_start'] - int(neighbor_data['Neighbor1_End'].iloc[0])), abs(row['subject_end'] - int(neighbor_data['Neighbor1_Start'].iloc[0])))
            if distance <= max_distance:
                blast_outcome = 1
                if row['subject_scaffold'] == neighbor_data['Neighbor2_Scaffold'].iloc[0] and neighbors_expected == 2:
                    distance2 = min(abs(row['subject_start'] - int(neighbor_data['Neighbor2_End'].iloc[0])), abs(row['subject_end'] - int(neighbor_data['Neighbor2_Start'].iloc[0])))
                    if distance2 <= max_distance:
                        blast_outcome = 2
                        # if both neighbors are present, we can stop searching
                        return(neighbors_expected, blast_outcome)
        blast_results.append(blast_outcome)
    return(neighbors_expected, max(blast_results))

def read_blast(file_path):
    # parse the BLAST output by skipping the first five lines, treating the rest as tab-separated values
    blastdf = pd.read_csv(file_path, sep='\t', comment='#', header=None)
    # rename columns
    blastdf.columns = [
        'query_seqid', 'subject_seqid', 'percent_identity', 'alignment_length', 'mismatches_count', 'gapopen_count',
        'query_start', 'query_end', 'subject_start', 'subject_end', 'evalue', 'bitscore', 'query_length', 'subject_length', 'subject_strand'
    ]
    return blastdf

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a BLAST output file in tabular format.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide a file name for the output file.''',
        default=None
        )
    parser.add_argument(
        '--neighbor_db','-ndb',type=str,
        help='''Provide a path to a directory containing neighbor data for accessory genes. Files should be in the format [isolate_name]_accessory_neighborfile.tsv.''',
        default=None
        )
    parser.add_argument(
        '--minimum_identity','-mi',type=float,
        help='''Provide the minimum percent identity for BLAST hits.''',
        default=0.9
        )
    parser.add_argument(
        '--minimum_evalue','-me',type=float,
        help='''Provide the minimum e-value for BLAST hits.''',
        default=1e-5
        )
    parser.add_argument(
        '--minimum_coverage','-mc',type=float,
        help='''Provide the minimum coverage for BLAST hits.''',
        default=0.9
        )
    parser.add_argument(
        '--maximum_distance','-md',type=float,
        help='''Provide the maximum distance from neighbors when evaluating BLAST hits.''',
        default=50000
        )
    args = parser.parse_args()
    blastdf = read_blast(args.input)
    outdata = evaluate_blast_hits_v2(
        blastdf, args.neighbor_db, minimum_identity=args.minimum_identity, minimum_evalue=args.minimum_evalue, 
        minimum_coverage=args.minimum_coverage, max_distance=args.maximum_distance)
    write_output_v2(outdata, args.output)

if __name__ == '__main__':
    main()
