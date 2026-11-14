
import argparse
import pandas as pd

def write_output_v2(outdata, output_file):
    with open(output_file, 'w') as out_f:
        _ = out_f.write('Query_SeqID\tSubject_SeqID\tExpected_Neighbors\tMax_Neighbors_in_BLAST\tGood_Blast_Hit\tOverlap_Found\tOverlap_Name\tOverlap_GeneFamily\tOverlap_Category\tPercent_Overlap\tPerfect_Match\tLong_Distance\n')
        for query_seqid, subject_hit_data_list in outdata.items():
            if subject_hit_data_list == []:
                print(f'Error: missing data for query {query_seqid}')
                quit(1)
            for subject_hit_data in subject_hit_data_list:
                out_line = query_seqid + '\t' + '\t'.join([str(x) for x in subject_hit_data]) + '\n'
                out_f.write(out_line)

def evaluate_blast_hits_v2(blastdf, neighbor_db_dir, minimum_identity=0.8, minimum_evalue=1e-5, minimum_coverage=0.8, 
    max_distance=50000, pangenome_matrix=None, filterfile=None, overlap_threshold=0.1, accessory_gene_file=None):
    # for each query sequence, return the list of subject sequence IDs that match the provided criteria
    # take all hits for a specific subject sequence ID into account, as long as they are above the thresholds
    pan = read_panaroo_data(pangenome_matrix, filterfile)
    isolate_list = [col for col in pan.columns if col != 'Gene']
    accessory_gene_df = pd.read_csv(accessory_gene_file)
    genes_not_seen = set(accessory_gene_df.iloc[:,0].tolist())
    outdata = {}
    for query_seqid, query_blastdf in blastdf.groupby('query_seqid'):
        if query_seqid not in genes_not_seen:
            print(f'warning, {query_seqid} has an unexpected name')
        genes_not_seen.discard(query_seqid)
        #print(f'Evaluating hits for query {query_seqid}')
        subject_hit_data_list = find_neighboring_best_hit_v2(
            query_blastdf, minimum_identity, minimum_evalue, minimum_coverage, max_distance, neighbor_db_dir, 
            isolate_list, overlap_threshold)
        outdata[query_seqid] = subject_hit_data_list
    # add each gene family that wasn't seen in the blast results
    for gf in genes_not_seen:
        subject_hit_data_list = []
        for isolate in isolate_list:
            subject_hit_data_list.append((isolate, 0, 0, False, False, 'None', 'None', 'None', 0, False, False))
        outdata[gf] = subject_hit_data_list
    return outdata

def find_neighboring_best_hit_v2(query_blastdf, minimum_identity, minimum_evalue, minimum_coverage, max_distance, neighbor_db_dir, isolate_list, overlap_threshold):
    # use a dictionary of neighboring gene families to identify if there is a high-quality hit in the expected location
    subject_hit_data = []
    # extract all hits that meet the criteria for each assembly
    #query_blastdf['assembly_name'] = query_blastdf['subject_seqid'].apply(lambda x: x.split('.scaffolds')[0])
    #query_blastdf['subject_scaffold'] = query_blastdf['subject_seqid'].apply(lambda x: x.split('.scaffolds_')[1])
    isolates_not_seen = set(isolate_list)
    for assembly_name, assembly_blastdf in query_blastdf.groupby('assembly_name'):
        if assembly_name not in isolate_list:
            continue
        query_gene_fam = query_blastdf['query_seqid'].iloc[0]
        assembly_blastdf['coverage'] = assembly_blastdf['alignment_length'] / assembly_blastdf['query_length']
        blastdf_filtered = assembly_blastdf[
            (assembly_blastdf['percent_identity'] >= minimum_identity * 100) &
            (assembly_blastdf['evalue'] <= minimum_evalue) &
            (assembly_blastdf['coverage'] >= minimum_coverage)
        ]
        neighbor_df = pd.read_csv(f'{neighbor_db_dir}/{assembly_name}_accessory_neighborfile.tsv', sep='\t')
        #neighbors_expected, max_neighbors_of_blast_hit, good_blast_hit = gff_neighbor_search_v2(query_gene_fam, assembly_name, blastdf_filtered, neighbor_df, max_distance)
        #subject_hit_data.append((assembly_name, neighbors_expected, max_neighbors_of_blast_hit, good_blast_hit))
        blast_data_summary = gff_neighbor_search_v4(query_gene_fam, assembly_name, blastdf_filtered, neighbor_df, max_distance, overlap_threshold)
        subject_hit_data.append(blast_data_summary)
        isolates_not_seen.discard(assembly_name)
    # for any isolates not seen at all in the blast results, add an entry indicating no hits
    for isolate in isolates_not_seen:
        subject_hit_data.append((isolate, 0, 0, False, False, 'None', 'None', 'None', 0, False, False))
    return subject_hit_data





def return_better_hit(current_hit_data, best_hit_data):
    if best_hit_data == None:
        return current_hit_data
    # assume that all blast hits are high-quality
    # if the current hit has more neighbors, always use it
    if current_hit_data['blast_neighbors'] > best_hit_data['blast_neighbors']:
        return current_hit_data
    # if the neighbors are the same, prefer the a hit that is a perfect match, followed by a hit that is not long distance
    elif current_hit_data['blast_neighbors'] == best_hit_data['blast_neighbors']:
        if best_hit_data['perfect_match'] == True:
            return best_hit_data
        if best_hit_data['perfect_match'] == False and current_hit_data['perfect_match'] == True:
            return current_hit_data
        if current_hit_data['long_distance'] == False and best_hit_data['long_distance'] == True:
            return current_hit_data
        # if long_distance is the same, prefer a hit without an overlap
        if current_hit_data['long_distance'] == best_hit_data['long_distance']:
            if current_hit_data['found_overlap'] == False and best_hit_data['found_overlap'] == True:
                return current_hit_data
    return best_hit_data



def gff_neighbor_search_v4(query_gene_fam, assembly_name, blastdf, neighbor_df, max_distance, overlap_threshold):
    # initialize outdata
    # isolate name, neighbors expected, max neighbors of best blast hit, if a good blast hit was found, if a good overlap was found, name of overlapping gene, name of overlapping gene family, category of overlapping gene
    outdata = (assembly_name, 0, 0, False, False, 'None', 'None', 'None', 0, False, False)
    # if no blast hits were found, return immediately
    if blastdf.empty:
        return outdata
    ## determine expected neighbors
    # subset df and ensure correct format
    neighbor_data = neighbor_df[neighbor_df['GeneFamily'] == query_gene_fam]
    if neighbor_data.empty:
        print(f'Error: no neighbor data found for gene family {query_gene_fam} in assembly {assembly_name}')
        quit(1)
    if neighbor_data['Neighbor1_Scaffold'].iloc[0] == 'absent' and neighbor_data['Neighbor2_Scaffold'].iloc[0] != 'absent':
        print('unexpected neighbor data format, neighbor 1 is absent but neighbor 2 is present!')
        quit(1)
    # count neighbors
    if neighbor_data['Neighbor1_Scaffold'].iloc[0] == 'absent':
        neighbors_expected = 0
    else:
        neighbors_expected = 1
        if neighbor_data['Neighbor2_Scaffold'].iloc[0] != 'absent':
            neighbors_expected = 2
    ## find best blast hit
    best_hit_data = None

    for _, row in blastdf.iterrows():
        # all rows in this table should consist of hits that pass the thresholds above
        # iterate through blast hits, recording the one with the most neighbors present
        # stop early if the number of neighbors is equal to the number of expected neighbors
        current_hit_data = {
        'blast_neighbors':0,'long_distance':False,'found_overlap':row['overlaps_annotation'],'overlap_gene':row['annotation_name'],
        'overlap_gene_family':row['annotation_gf'],'overlap_category':row['annotation_type'],'percent_overlap':round(row['annotation_coverage']/row['alignment_length'], 2),
        'perfect_match':False,'percent_identity':row['percent_identity'],'alignment_coverage':row['alignment_length']/row['query_length']
        }

        # check the first neighbor
        if row['subject_scaffold'] == neighbor_data['Neighbor1_Scaffold'].iloc[0] and neighbors_expected >= 1:
            current_hit_data['blast_neighbors'] += 1
            # check distance
            distance = min(abs(row['subject_start'] - int(neighbor_data['Neighbor1_End'].iloc[0])), abs(row['subject_end'] - int(neighbor_data['Neighbor1_Start'].iloc[0])))
            if distance >= max_distance:
                current_hit_data['long_distance'] = True
        # check the second neighbor
        if row['subject_scaffold'] == neighbor_data['Neighbor2_Scaffold'].iloc[0] and neighbors_expected == 2:
            current_hit_data['blast_neighbors'] += 1
            # check distance
            distance2 = min(abs(row['subject_start'] - int(neighbor_data['Neighbor2_End'].iloc[0])), abs(row['subject_end'] - int(neighbor_data['Neighbor2_Start'].iloc[0])))
            if distance2 >= max_distance:
                current_hit_data['long_distance'] = True

        # determine if this is a perfect match - 100% identity and 100% coverage of the query
        if current_hit_data['percent_identity'] == 100 and current_hit_data['alignment_coverage'] == 1:
            current_hit_data['perfect_match'] = True

        best_hit_data = return_better_hit(current_hit_data, best_hit_data)

    # if the best hit has an overlap below the threshold, remove it
    if best_hit_data['percent_overlap'] < overlap_threshold and best_hit_data['found_overlap'] == True:
        best_hit_data['found_overlap'] = False
        best_hit_data['overlap_gene'] = 'None'
        best_hit_data['overlap_gene_family'] = 'None'
        best_hit_data['overlap_category'] = 'None'
        best_hit_data['percent_overlap'] = 0

    # return data for the best hit
    outdata = (
        assembly_name, neighbors_expected, best_hit_data['blast_neighbors'], True, 
        best_hit_data['found_overlap'], best_hit_data['overlap_gene'], best_hit_data['overlap_gene_family'], best_hit_data['overlap_category'],
        best_hit_data['percent_overlap'], best_hit_data['perfect_match'], best_hit_data['long_distance']
        )
    return outdata






def read_panaroo_data(panaroo_file, filterfile):
    pan = pd.read_csv(panaroo_file,keep_default_na=False)
    # subset pan to only contain columns for the gene family names and isolates
    pan = pan.loc[:, pan.columns.str.startswith(("Gene", "SRR", "ARR", "DRR", "ERR", "UM_", "Chi_", "b8441"))]
    # subset pan to only include isolates in the filter file, if provided
    if filterfile is not None:
        with open(filterfile, 'r') as fh:
            filter_list = [line.strip() for line in fh]
        # always keep the 'Gene' column
        filter_list.append('Gene')
        # include only the columns in the filter list
        pan = pan.loc[:, pan.columns.isin(filter_list)]
    return pan

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
        '--maximum_evalue','-me',type=float,
        help='''Provide the maximum e-value for BLAST hits.''',
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
    parser.add_argument(
        '--pangenome','-pan',type=str,
        help='''Provide the path to a pangenome matrix. This should contain names for genes.''',
        default=None
        )
    parser.add_argument(
        '--filterfile','-ff',type=str,
        help='''Provide the path to the filter file containing isolate names. The pangenome matrix will be subset to only include these isolates.''',
        default=None
        )
    parser.add_argument(
        '--overlap_threshold','-ot',type=float,
        help='''Provide the minimum overlap threshold for evaluating BLAST hits.''',
        default=0.1
        )
    parser.add_argument(
        '--gene_list','-gl',type=str,
        help='''Provide the path to a one-column csv containing a list of accessory genes.''',
        default=None
        )
    args = parser.parse_args()
    blastdf = pd.read_csv(args.input, sep='\t', keep_default_na=False)
    outdata = evaluate_blast_hits_v2(
        blastdf, args.neighbor_db, minimum_identity=args.minimum_identity, minimum_evalue=args.maximum_evalue, 
        minimum_coverage=args.minimum_coverage, max_distance=args.maximum_distance, pangenome_matrix=args.pangenome, filterfile=args.filterfile, 
        overlap_threshold=args.overlap_threshold, accessory_gene_file=args.gene_list)
    write_output_v2(outdata, args.output)

if __name__ == '__main__':
    main()
