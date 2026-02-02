
import argparse
import pandas as pd
from itertools import filterfalse

def read_blast(file_path):
    # parse the BLAST output by skipping the first five lines, treating the rest as tab-separated values
    blastdf = pd.read_csv(file_path, sep='\t', comment='#', header=None)
    # rename columns
    blastdf.columns = [
        'query_seqid', 'subject_seqid', 'percent_identity', 'alignment_length', 'mismatches_count', 'gapopen_count',
        'query_start', 'query_end', 'subject_start', 'subject_end', 'evalue', 'bitscore', 'query_length', 'subject_length', 'subject_strand'
    ]
    return blastdf

def find_best_hit(query_blastdf, minimum_identity, minimum_evalue, minimum_coverage):
    #print(f'processing {query_blastdf["query_seqid"].iloc[0]}')
    subject_hits = []
    # extract the single best hit that meets the criteria for each assembly
    query_blastdf['assembly_name'] = query_blastdf['subject_seqid'].apply(lambda x: x.split('.scaffolds')[0])
    for assembly_name, assembly_blastdf in query_blastdf.groupby('assembly_name'):
        #print(f'checking assembly {assembly_name}')
        assembly_blastdf['coverage'] = assembly_blastdf['alignment_length'] / assembly_blastdf['query_length']
        blastdf_filtered = assembly_blastdf[
            (assembly_blastdf['percent_identity'] >= minimum_identity * 100) &
            (assembly_blastdf['evalue'] <= minimum_evalue) &
            (assembly_blastdf['coverage'] >= minimum_coverage)
        ]
        if blastdf_filtered.empty:
            continue
        #print(blastdf_filtered)
        # find the hit with the highest coverage
        best_hit = blastdf_filtered.loc[blastdf_filtered['coverage'].idxmax()]
        #print(best_hit)
        assembly_name = best_hit['assembly_name']
        if assembly_name not in subject_hits:
            subject_hits.append(assembly_name)
    #print(subject_hits)
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


def evaluate_blast_hits(blastdf, minimum_identity=0.3, minimum_evalue=1e-5, minimum_coverage=0.8, summary_strategy='best_hit'):
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
        outdata[query_seqid] = subject_hits
    return outdata

def write_output(outdata, output_file):
    with open(output_file, 'w') as out_f:
        _ = out_f.write('Query_SeqID\tSubject_SeqID\n')
        for query_seqid, subject_hits in outdata.items():
            if subject_hits == []:
                subject_hits = ['None']
            for subject_hit in subject_hits:
                out_f.write(f'{query_seqid}\t{subject_hit}\n')

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
        '--strategy','-s',type=str,
        help='''Specify a strategy for summarizing BLAST hits. best_hit will take only the best hit from each subject.
        merge_hits will merge all hits from the same chromosome together before calculating coverage.''',
        default='best_hit'
        )    
    args = parser.parse_args()
    blastdf = read_blast(args.input)
    outdata = evaluate_blast_hits(blastdf)
    #print(outdata)
    write_output(outdata, args.output)

if __name__ == '__main__':
    main()
