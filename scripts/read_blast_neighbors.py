
import argparse
import pandas as pd
from itertools import filterfalse
import networkx as nx
import gffutils as gff

# current plan:
# use the pangenome graph to identify the nearest upstream and downstream genes
    # this will be different depending on the isolate
    # for now, only consider cases where both of the nearest flanking gene families are on the same contig
# perform an additional search for BLAST hits, subsetting to the hits that are close to this genomic region


def read_blast(file_path):
    # parse the BLAST output by skipping the first five lines, treating the rest as tab-separated values
    blastdf = pd.read_csv(file_path, sep='\t', comment='#', header=None)
    # rename columns
    blastdf.columns = [
        'query_seqid', 'subject_seqid', 'percent_identity', 'alignment_length', 'mismatches_count', 'gapopen_count',
        'query_start', 'query_end', 'subject_start', 'subject_end', 'evalue', 'bitscore', 'query_length', 'subject_length', 'subject_strand'
    ]
    return blastdf

def read_pangenome_graph(graph_path):
    pangraph = nx.read_gml(graph_path, label ='name')
    neighbor_dict = {}
    for node in pangraph.nodes:
        neighbors = list(pangraph.neighbors(node))
        neighbor_dict[node] = neighbors
    return neighbor_dict

def read_panaroo_data(panaroo_file, filterfile=None):
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

def find_neighboring_best_hit(query_blastdf, minimum_identity, minimum_evalue, minimum_coverage, max_distance, neighbor_dict, gff_dir, pan_matrix):
    # use a dictionary of neighboring gene families to identify if there is a high-quality hit in the expected location
    neighbor_present_hits, neighbor_absent_hits, no_neighbor_hits = [], [], []
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
        # add a column that indicates if the hit is within max_distance of the expected region
        neighbor_status = gff_neighbor_search(blastdf_filtered, neighbor_dict, gff_dir, max_distance, pan_matrix)
        if neighbor_status == 'present' and assembly_name not in neighbor_present_hits:
            neighbor_present_hits.append(assembly_name)
        elif neighbor_status == 'absent' and assembly_name not in neighbor_absent_hits:
            neighbor_absent_hits.append(assembly_name)
        elif neighbor_status == 'no_neighbors' and assembly_name not in no_neighbor_hits:
            no_neighbor_hits.append(assembly_name)
    return neighbor_present_hits, neighbor_absent_hits, no_neighbor_hits

def gff_neighbor_search(blastdf, neighbor_dict, gff_dir, max_distance, pan_matrix):
    # blastdf should consist of a single gene family query and a single assembly name
    # first, extract the gene family and the assembly name
    query_gene_fam = blastdf['query_seqid'].iloc[0]
    subject_assembly = blastdf['assembly_name'].iloc[0]
    #print(f'Checking neighbors for gene family {query_gene_fam} in assembly {subject_assembly}')
    # get the neighboring gene families for this query gene family
    neighbor_list = neighbor_dict.get(query_gene_fam, [])
    #print(f'initial neighbor list: {neighbor_list}')
    present_neighbor_list = []
    for neighbor in neighbor_list:
        # check if these neighbors are actually present in this assembly, based on the pangenome matrix
        # if they are not, use the next neighbor (if one exists)
        target_genefam = neighbor
        genes_seen = []
        # b8441 coding sequences start with cds- and end in .1
        # all other isolate coding sequences end in .cds
        while not pan_matrix.loc[pan_matrix['Gene'] == target_genefam, subject_assembly].iloc[0].endswith(('.cds','.1')):
            # note: the above condition does not currently handle cases where there is a list of multiple genes that ends with a refound or pseudo gene
            # example: 'FUN_001884-T1.cds;FUN_001883-T1.cds_pseudo'
            genes_seen.append(target_genefam)
            next_neighbors = neighbor_dict.get(target_genefam, [])
            # remove any previously seen genes to avoid infinite loops
            next_neighbors = [x for x in next_neighbors if x not in genes_seen]
            if next_neighbors == [] or len(genes_seen) > 100:
                target_genefam = None
                break
            # if there are multiple next neighbors, just take the first one
            target_genefam = next_neighbors[0]
        if target_genefam is not None:
            cds_string = pan_matrix.loc[pan_matrix['Gene'] == target_genefam, subject_assembly].iloc[0]
            # this can be multiple genes, separated by a semicolon
            cds_list = [x for x in cds_string.split(';') if x.endswith(('.cds','.1'))]
            if cds_list != []:
                present_neighbor_list.append(cds_list[0])
    # the neighbor list should now consist only of neighboring gene families that are actually present in this assembly
    #print(f'present neighbor list: {present_neighbor_list}')
    # for each of these neighbors, extract their scaffold, start, and end positions from the gff file
    # note: this is probably too inefficient for a large number of queries, since this will recreate the gffutils db each time
    # consider making a single large gffutils db for all assemblies instead in the future
    if present_neighbor_list == []:
        return 'no_neighbors'
    else:
        neighbor_positions = []
        gff_file = f'{gff_dir}/{subject_assembly}.gff'
        gff_db = gff.create_db(gff_file,dbfn=":memory:",force=True,keep_order=False,merge_strategy="create_unique",sort_attribute_values=True,from_string=False)
        for neighbor_gene in present_neighbor_list:
            gene_feature = gff_db[neighbor_gene]
            scaffold = gene_feature.seqid
            start = gene_feature.start
            end = gene_feature.end
            neighbor_positions.append((scaffold, start, end))
    #print(f'neighbor positions: {neighbor_positions}')
    # for each neighbor position, check if it is on the same scaffold and within the max_distance of any of the BLAST hits
    # for now, return a hit if any BLAST hit is in close proximity to any neighbor
    for _, row in blastdf.iterrows():
        scaffold1, start1, end1 = neighbor_positions[0]
        if row['subject_scaffold'] == scaffold1:
            # check distance
            distance = min(abs(row['subject_start'] - end1), abs(row['subject_end'] - start1))
            if distance <= max_distance:
                if len(neighbor_positions) == 1:
                    #print(f'found valid BLAST hit near neighbor: {row}')
                    return 'present'
                else: # check the second neighbor as well
                    scaffold2, start2, end2 = neighbor_positions[1]
                    distance2 = min(abs(row['subject_start'] - end2), abs(row['subject_end'] - start2))
                    if distance2 <= max_distance:
                        #print(f'found valid BLAST hit near both neighbors: {row}')
                        return 'present'
    return 'absent'
                

def find_best_hit(query_blastdf, minimum_identity, minimum_evalue, minimum_coverage):
    subject_hits = []
    # extract the single best hit that meets the criteria for each assembly
    query_blastdf['assembly_name'] = query_blastdf['subject_seqid'].apply(lambda x: x.split('.scaffolds')[0])
    for assembly_name, assembly_blastdf in query_blastdf.groupby('assembly_name'):
        assembly_blastdf['coverage'] = assembly_blastdf['alignment_length'] / assembly_blastdf['query_length']
        blastdf_filtered = assembly_blastdf[
            (assembly_blastdf['percent_identity'] >= minimum_identity * 100) &
            (assembly_blastdf['evalue'] <= minimum_evalue) &
            (assembly_blastdf['coverage'] >= minimum_coverage)
        ]
        if blastdf_filtered.empty:
            continue
        # find the hit with the highest coverage, if needed
        #best_hit = blastdf_filtered.loc[blastdf_filtered['coverage'].idxmax()]
        #print(best_hit)
        if assembly_name not in subject_hits:
            subject_hits.append(assembly_name)
    return subject_hits

def find_merged_hits(query_blastdf, minimum_identity, minimum_evalue, minimum_coverage):
    # Important Note: this does not consider if the hits are continuous or adjacent, so fragmented hits spread out across 
    # a single chromosome will still count towards coverage
    subject_hits = []
    query_length = query_blastdf['query_length'].iloc[0]
    for subject_seqid, subject_blastdf in query_blastdf.groupby('subject_seqid'):
        # subject_seqid contains both the assembly name and the contig/chromosome name
        blastdf_filtered = subject_blastdf[
            (subject_blastdf['percent_identity'] >= minimum_identity * 100) &
            (subject_blastdf['evalue'] <= minimum_evalue)
        ]
        if blastdf_filtered.empty:
            continue
        #print(f'Checking hits for subject {subject_seqid} against query {query_seqid}')
        #print(blastdf_filtered)
        # calculate how much of the query sequence is covered by these hits
        # do this by subtracting each range from the full length of the query sequence
        remaining_query = range(1, query_length + 1)
        for _, row in blastdf_filtered.iterrows():
            match_range = range(row['query_start'], row['query_end'] + 1)
            remaining_query = filterfalse(lambda x: x in match_range, remaining_query)
        # any remaining positions should not be present in any of the hits
        non_covered_length = len(list(remaining_query))
        coverage = (query_length - non_covered_length) / query_length
        #print(f'Calculated coverage of {coverage:.2f}')
        if coverage >= minimum_coverage:
            #print(f'Subject passed')
            # eventually, I'd like to store information about the quality of the hits here
            # for now, just save the assembly name, removing the contig/chromosome name
            assembly_name = subject_seqid.split('.scaffolds')[0]
            if assembly_name not in subject_hits:
                subject_hits.append(assembly_name)
    return subject_hits


def evaluate_blast_hits(blastdf, minimum_identity=0.3, minimum_evalue=1e-5, minimum_coverage=0.8, summary_strategy='best_hit', 
    max_distance=5000, neighbor_dict={}, gff_dir='', pan_matrix=pd.DataFrame()):
    # for each query sequence, return the list of subject sequence IDs that match the provided criteria
    # take all hits for a specific subject sequence ID into account, as long as they are above the thresholds
    outdata = {}
    for query_seqid, query_blastdf in blastdf.groupby('query_seqid'):
        #print(f'Evaluating hits for query {query_seqid}')
        subject_hits = []
        if summary_strategy == 'best_hit':
            subject_hits = find_best_hit(query_blastdf, minimum_identity, minimum_evalue, minimum_coverage)
        if summary_strategy == 'merge_hits':
            subject_hits = find_merged_hits(query_blastdf, minimum_identity, minimum_evalue, minimum_coverage)
        if summary_strategy == 'neighbor_best_hit':
            neighbor_present_hits, neighbor_absent_hits, no_neighbor_hits = find_neighboring_best_hit(query_blastdf, minimum_identity, minimum_evalue, minimum_coverage, max_distance, neighbor_dict, gff_dir, pan_matrix)
            subject_hits = [neighbor_present_hits, neighbor_absent_hits, no_neighbor_hits]
        outdata[query_seqid] = subject_hits
    return outdata

def write_output(outdata, output_file, summary_strategy='best_hit'):
    with open(output_file, 'w') as out_f:
        if summary_strategy == 'best_hit' or summary_strategy == 'merge_hits':
            _ = out_f.write('Query_SeqID\tSubject_SeqID\n')
            for query_seqid, subject_hits in outdata.items():
                if subject_hits == []:
                    subject_hits = ['None']
                for subject_hit in subject_hits:
                    out_f.write(f'{query_seqid}\t{subject_hit}\n')
        if summary_strategy == 'neighbor_best_hit':
            _ = out_f.write('Query_SeqID\tSubject_SeqID\tSubject_Neighbor_Status\n')
            for query_seqid, subject_hits in outdata.items():
                if subject_hits == [[], [], []]:
                    subject_hits = [[], [], ['None']]
                for neighbor_present_hit in subject_hits[0]:
                    out_f.write(f'{query_seqid}\t{neighbor_present_hit}\tpresent\n')
                for neighbor_absent_hit in subject_hits[1]:
                    out_f.write(f'{query_seqid}\t{neighbor_absent_hit}\tabsent\n')
                for no_neighbor_hit in subject_hits[2]:
                    out_f.write(f'{query_seqid}\t{no_neighbor_hit}\tno_neighbors\n')

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
        '--strategy','-s',choices=['best_hit','merge_hits','neighbor_best_hit'],type=str,
        help='''Specify a strategy for summarizing BLAST hits. best_hit will take only the best hit from each subject.
        merge_hits will merge all hits from the same chromosome together before calculating coverage.''',
        default='best_hit'
        )
    parser.add_argument(
        '--graph','-gml',type=str,
        help='''Provide a path to a pangenome graph file in GML format. Only needed for the neighbor_best_hit strategy.''',
        default=None
        )
    parser.add_argument(
        '--matrix','-m',type=str,
        help='''Provide a path to a pangenome matrix file, containing CDS names. Only needed for the neighbor_best_hit strategy.''',
        default=None
        )
    parser.add_argument(
        '--gffdir','-gff',type=str,
        help='''Provide a path to a directory of GFF files. These should match the isolate names in the pangenome data. Only needed for the neighbor_best_hit strategy.''',
        default=None
        )
    args = parser.parse_args()
    blastdf = read_blast(args.input)
    if args.strategy == 'neighbor_best_hit':
        neighbor_dict = read_pangenome_graph(args.graph)
        pan_matrix = read_panaroo_data(args.matrix)
        outdata = evaluate_blast_hits(
            blastdf, summary_strategy=args.strategy, neighbor_dict=neighbor_dict,
            gff_dir=args.gffdir, pan_matrix=pan_matrix
        )
    else:
        outdata = evaluate_blast_hits(blastdf, args.strategy)
    #print(outdata)
    write_output(outdata, args.output, args.strategy)

if __name__ == '__main__':
    main()
