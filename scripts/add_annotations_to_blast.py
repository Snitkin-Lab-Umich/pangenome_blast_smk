
import argparse
import pandas as pd
import gffutils as gff
import os

def read_blast(file_path):
    # parse the BLAST output by skipping the first five lines, treating the rest as tab-separated values
    blastdf = pd.read_csv(file_path, sep='\t', comment='#', header=None)
    # rename columns
    blastdf.columns = [
        'query_seqid', 'subject_seqid', 'percent_identity', 'alignment_length', 'mismatches_count', 'gapopen_count',
        'query_start', 'query_end', 'subject_start', 'subject_end', 'evalue', 'bitscore', 'query_length', 'subject_length', 'subject_strand'
    ]
    return blastdf

def read_panaroo_data(panaroo_file, filterfile, keep_cols=False):
    pan = pd.read_csv(panaroo_file,keep_default_na=False)
    # subset pan to only contain columns for the gene family names and isolates
    if not keep_cols:
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

def interval_overlap(s1,e1,s2,e2):
    if s1 > e1 or s2 > e2:
        print('interval input error')
        quit(1)
    ovs = max(s1,s2)
    ove = min(e1,e2)
    ov = max(0,(ove-ovs+1))
    return(ov)

def get_overlapping_annotations(gffdb, cds_to_gf_dict, gf_to_category_dict, assembly_name, subject_scaffold, subject_start, subject_end):
    # find any CDS annotations in gffdb that overlap the provided interval
    # if there are multiple, choose the one with the longest overlap
    # return the gene family of the chosen annotation, if it exists
    best_overlap = (False, 0, 'None', 'None', 'None') # initialize as no overlap, 0 overlap length, and no annotation info
    # search through all CDS and mRNA annotations
    # for ann in gffdb.region(region=(subject_scaffold, subject_start, subject_end), featuretype=['CDS', 'mRNA'], completely_within=False):
    # search through only CDS annotations
    for ann in gffdb.region(region=(subject_scaffold, subject_start, subject_end), featuretype=['CDS'], completely_within=False):
        overlap = interval_overlap(subject_start, subject_end, ann.start, ann.end)
        if overlap > best_overlap[1]:
            # find this annotation's cds name in pan_name to get the gene family
            # if this annotation is an mRNA, find its first child CDS annotation and search for that in pan_name
            if ann.featuretype == 'CDS':
                cds_name = ann.id
            elif ann.featuretype == 'mRNA':
                children = list(gffdb.children(ann.id, featuretype='CDS', order_by='start'))
                if len(children) > 0:
                    cds_name = children[0].id
                else:
                    cds_name = ann.id
            if cds_name in cds_to_gf_dict[assembly_name]:
                gene_family = cds_to_gf_dict[assembly_name][cds_name]
                gene_category = gf_to_category_dict[gene_family]
                best_overlap = (True, overlap, cds_name, gene_family, gene_category)
            else:
                best_overlap = (True, overlap, cds_name, 'None', 'None')
    return best_overlap

def make_cds_to_gf_dict(pan_name):
    pan = pan_name.loc[:, pan_name.columns.str.startswith(("Gene", "SRR", "ARR", "DRR", "ERR", "UM_", "Chi_", "b8441"))]
    cds_to_gf_dict = {}
    for isolate in pan_name.columns:
        if isolate == 'Gene':
            continue
        cds_to_gf_dict[isolate] = {}
    for index, row in pan.iterrows():
        gene_family = pan.loc[index, 'Gene']
        for isolate in pan.columns:
            if isolate == 'Gene':
                continue
            cds_string = pan.loc[index, isolate].split(',')
            for cds in cds_string:
                if cds != '' and 'refound' not in cds and 'pseudo' not in cds:
                    cds_to_gf_dict[isolate][cds] = gene_family
    return cds_to_gf_dict

def make_gf_to_category_dict(pan_num, clade_col):
    pan = pan_num.loc[:, pan_num.columns.str.startswith(("Gene", clade_col))]
    gf_to_category_dict = {}
    for index, row in pan.iterrows():
        gene_family = pan.loc[index, 'Gene']
        category = pan.loc[index, clade_col]
        gf_to_category_dict[gene_family] = category
    return gf_to_category_dict

def add_annotations(blast_file, gff_dir, pangenome_name_matrix, pangenome_num_matrix, clade_name, output_file):
    if clade_name == 'clade_I':
        clade_col = 'category_cI'
    elif clade_name == 'clade_III':
        clade_col = 'category_cIII'
    elif clade_name == 'clade_IV':
        clade_col = 'category_cIV'
    else:
        print('Unrecognized clade name, defaulting to species-wide categories')
        clade_col = 'category'
    blastdf = read_blast(blast_file)
    # read in both pangenome dataframes
    pan_name = read_panaroo_data(pangenome_name_matrix, filterfile=None, keep_cols=False)
    pan_num = read_panaroo_data(pangenome_num_matrix, filterfile=None, keep_cols=True)
    # make dictionaries for both
    cds_to_gf_dict = make_cds_to_gf_dict(pan_name)
    gf_to_category_dict = make_gf_to_category_dict(pan_num, clade_col)
    # add columns to blastdf for annotation info
    blastdf['assembly_name'] = blastdf['subject_seqid'].apply(lambda x: x.split('.scaffolds')[0])
    blastdf['subject_scaffold'] = blastdf['subject_seqid'].apply(lambda x: x.split('.scaffolds_')[1])
    # iterate through each unique isolate in blastdf's assembly_name column
    # on each iteration, update the rows of blastdf with that assembly name to include annotation data
    blastdf['overlaps_annotation'] = False
    blastdf['annotation_coverage'] = 0
    blastdf['annotation_name'] = 'None'
    blastdf['annotation_gf'] = 'None'
    blastdf['annotation_type'] = 'None'
    for assembly_name in blastdf['assembly_name'].unique():
        gff_file_path = os.path.join(gff_dir, f'{assembly_name}.gff')
        if not os.path.isfile(gff_file_path):
            gff_file_path = os.path.join(gff_dir, f'{assembly_name}.gff3')
        if not os.path.isfile(gff_file_path):
            print(f'Unable to locate gff file for assembly {assembly_name} at {gff_file_path}')
            continue
        gffdb = gff.create_db(gff_file_path, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
        for index, row in blastdf[blastdf['assembly_name'] == assembly_name].iterrows():
            subject_scaffold = row['subject_scaffold']
            subject_blast_start = int(row['subject_start'])
            subject_blast_end = int(row['subject_end'])
            subject_start = min(subject_blast_start,subject_blast_end)
            subject_end = max(subject_blast_start,subject_blast_end)
            overlap_data = get_overlapping_annotations(gffdb, cds_to_gf_dict, gf_to_category_dict, assembly_name, subject_scaffold, subject_start, subject_end)
            blastdf.at[index, 'overlaps_annotation'] = overlap_data[0]
            blastdf.at[index, 'annotation_coverage'] = overlap_data[1]
            blastdf.at[index, 'annotation_name'] = overlap_data[2]
            blastdf.at[index, 'annotation_gf'] = overlap_data[3]
            blastdf.at[index, 'annotation_type'] = overlap_data[4]
    blastdf.to_csv(output_file, sep='\t', index=False)



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
        '--gff_dir','-gff',type=str,
        help='''Provide a path to a directory containing GFF files for each isolate. Files should be in the format [isolate_name].gff or [isolate_name].gff3.''',
        default=None
        )
    parser.add_argument(
        '--pangenome_names','-pname',type=str,
        help='''Provide the path to a pangenome matrix. This should contain names for genes.''',
        default=None
        )
    parser.add_argument(
        '--pangenome_numbers','-pnum',type=str,
        help='''Provide the path to a pangenome matrix. This should contain numbers for genes.''',
        default=None
        )
    parser.add_argument(
        '--clade_name','-c',type=str,
        help='''Provide the name of the clade. This should correspond to a category column in the pangenome matrix.''',
        default=None
        )
    args = parser.parse_args()
    add_annotations(args.input, args.gff_dir, args.pangenome_names, args.pangenome_numbers, args.clade_name, args.output)

if __name__ == '__main__':
    main()
