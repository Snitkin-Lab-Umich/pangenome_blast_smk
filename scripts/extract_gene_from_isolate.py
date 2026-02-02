import os
import subprocess
from Bio import SeqIO
import argparse
import gffutils as gff
import pandas as pd

def extract_sequences(input_file, assemblies_dir, gff_dir, pangenome_file, filterfile, output_file):
    # steps for each gene family
    # 1) read the pangenome matrix, get the row that matches the gene family name
    # 2) get the list of isolates that actually have that gene family
    # 3) for each isolate, parse the gff file to get the start and end positions of that gene
    # 4) extract this sequence from the corresponding fasta file
    # 5) compare the extracted sequences and write the longest/best one to the output file
    #####
    genefam_data = pd.read_csv(input_file)
    fullpan = pd.read_csv(pangenome_file,keep_default_na=False)
    # subset pan to only contain columns for the gene family names and isolates
    fullpan = fullpan.loc[:, fullpan.columns.str.startswith(("Gene", "SRR", "ARR", "DRR", "ERR", "UM_", "Chi_", "b8441"))]
    # subset pan to only include isolates in the filter file, if provided
    if filterfile is not None:
        with open(filterfile, 'r') as fh:
            filter_list = [line.strip() for line in fh]
        # always keep the 'Gene' column
        filter_list.append('Gene')
        # include only the columns in the filter list
        pan = fullpan.loc[:, fullpan.columns.isin(filter_list)]
    with open(output_file, 'w') as out_fasta:
        for genefam in genefam_data.iloc[:,0]:
            # if genefam not in pan['Gene'].values:
            #     print(f'{genefam} not found in pangenome file')
            #     continue
            # genefam_row_index = pan[pan['Gene'] == genefam].index[0]
            # pan_data_series = pan.loc[genefam_row_index]
            # # create a dictionary of isolate:gene_name for the isolates that have this gene family
            # isolate_gene_dict = {}
            # for isolate in pan_data_series.index:
            #     if isolate == 'Gene':
            #         continue
            #     gene_name = pan_data_series[isolate]
            #     if gene_name == '' or 'refound' in gene_name or 'pseudo' in gene_name:
            #         continue
            #     isolate_gene_dict[isolate] = gene_name
            isolate_gene_dict = make_isolate_dict(genefam, pan)
            if isolate_gene_dict == {}:
                print(f'Warning: no searchable genes found in gene family {genefam} for this clade')
                # use the full pan data to recover the sequence, if it cannot be found within the clade
                isolate_gene_dict = make_isolate_dict(genefam, fullpan)
            if isolate_gene_dict == {}:
                print(f'Warning: no searchable genes found in gene family {genefam} at all')
                continue
            isolate, gene_name = find_best_gene(isolate_gene_dict, gff_dir)
            # remove the .cds suffix from the gene name to extract the mRNA feature rather than the CDS feature
            # if using the spliced data, use the CDS name instead
            #gene_name_mrna = gene_name.replace('.cds', '')
            gene_name_mrna = gene_name
            gff_file = os.path.join(gff_dir, f'{isolate}.gff3')
            if not os.path.isfile(gff_file):
                gff_file = os.path.join(gff_dir, f'{isolate}.gff')
            if not os.path.isfile(gff_file):
                print(f'Unable to locate gff file for isolate {isolate} at {gff_file}')
                continue
            fasta_file = os.path.join(assemblies_dir, f'{isolate}.scaffolds.fa')
            query_fasta_record = extract_record_from_gff(gff_file, fasta_file, gene_name_mrna)
            #print(f'Found gene {query_fasta_record.id} in isolate {isolate} for gene family {genefam}')
            #print(f'gff: {gff_file}, fasta: {fasta_file}, gene_name: {gene_name_mrna}')
            # write the record to the output fasta file
            query_fasta_record.id = genefam
            query_fasta_record.description = ''
            SeqIO.write(query_fasta_record, out_fasta, 'fasta')

def make_isolate_dict(genefam, pan):
    if genefam not in pan['Gene'].values:
        print(f'{genefam} not found in pangenome file')
        return {}
    genefam_row_index = pan[pan['Gene'] == genefam].index[0]
    pan_data_series = pan.loc[genefam_row_index]
    # create a dictionary of isolate:gene_name for the isolates that have this gene family
    isolate_gene_dict = {}
    for isolate in pan_data_series.index:
        if isolate == 'Gene':
            continue
        gene_name = pan_data_series[isolate]
        if ';' in gene_name:
            gene_name_list = gene_name.split(';')
        else:
            gene_name_list = [gene_name]
        # keep only genes that end with either .cds or .1
        gene_name_list = [x for x in gene_name_list if x.endswith(('.cds','.1'))]
        gene_name_final = ';'.join(gene_name_list)
        if gene_name_final == '':
            continue
        isolate_gene_dict[isolate] = gene_name_final
        # original method to explicitly exclude pseudo genes
        # if gene_name == '' or 'refound' in gene_name or 'pseudo' in gene_name:
        #     continue
        # method to keep pseudo genes in the search
        # if gene_name.endswith('_pseudo'):
        #     gene_name = gene_name.replace('_pseudo','')
        # isolate_gene_dict[isolate] = gene_name
    return isolate_gene_dict

def extract_record_from_gff(gff_file, fasta_file, gene_name):
    # parse the gff file to get the start and end positions of the gene
    db = gff.create_db(gff_file, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    gene_feature = db[gene_name]
    seqid = gene_feature.seqid
    start = gene_feature.start
    end = gene_feature.end
    strand = gene_feature.strand
    # extract the sequence from the fasta file
    #print(f'Extracting sequence for gene {gene_name} from {fasta_file}, seqid: {seqid}, start: {start}, end: {end}, strand: {strand}')
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id == seqid:
            gene_seq = record.seq[start-1:end]  # gff is 1-based, python is 0-based
            gene_record = SeqIO.SeqRecord(gene_seq, id=gene_name, description='')
            return(gene_record)
    print(f'Could not find sequence for gene {gene_name} in fasta file {fasta_file}')
    quit(1)



def find_best_gene(isolate_gene_dict, gff_dir):
    # take a dict of isolate:gene_name
    # eventually, this will need to compare lengths/quality of sequences
    # for now, simply prioritize any isolates that start with UM_ or Chi_, as these are hybrid assemblies
    final_isolate = None
    final_gene_name = None
    for isolate in isolate_gene_dict:
        if isolate.startswith('UM_') or isolate.startswith('Chi_'):
            final_isolate = isolate
            final_gene_name = isolate_gene_dict[isolate]
            break
        final_isolate = isolate
        final_gene_name = isolate_gene_dict[isolate]
    # some gene names are actually lists of multiple genes separated by a semicolon
    if ';' in final_gene_name:
        final_gene_name = final_gene_name.split(';')[0]
    return final_isolate, final_gene_name


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a one-column csv containing the names of accessory gene families to extract sequences for. These should be the exact group names from Panaroo.''',
        default=None
        )
    parser.add_argument(
        '--assemblies','-a',type=str,
        help='''Provide the name of a directory containing assemblies for all isolates. This should be in the same format as the Panaroo input, with
        file names ending in .scaffolds.fa''',
        default=None
        )
    parser.add_argument(
        '--gff','-g',type=str,
        help='''Provide the name of a directory containing gff files for all isolates. This should be in the same format as the Panaroo input, with
        file names ending in .gff3''',
        default=None
        )
    parser.add_argument(
        '--pangenome','-p',type=str,
        help='''Provide a path to the presence/absence matrix provided by Panaroo. Gene names should be present in this file.''',
        default=None
        )
    parser.add_argument(
        '--filterfile','-ff',type=str,
        help='''Provide a path to a one-column csv consisting only of isolate names. Only isolates in this file will be used.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide the path to the output file. This will be a fasta file with a single entry for each gene family in the input list.''',
        default=None
        )
    args = parser.parse_args()
    extract_sequences(args.input, args.assemblies, args.gff, args.pangenome, args.filterfile,args.output)


if __name__ == '__main__':
    main()