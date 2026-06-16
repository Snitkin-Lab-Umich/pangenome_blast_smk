
import argparse
import pandas as pd

# summarize the results of the blast comparison as best as possible
# compare the number of isolates with each type of blast hit to the isolates that are present or absent in the panaroo data

def compare_panaroo_blast_v3(panaroo_file, blast_file, filterfile, output_file):
    # read in the panaroo presence/absence matrix
    pandf = read_panaroo_data(panaroo_file, filterfile)
    # read in the blast results
    blastdf = pd.read_csv(blast_file, sep='\t', keep_default_na=False)
    # look at each unique query_seqid in blastdf and compare to panaroo
    with open(output_file, 'w') as fh:
        header_line = 'Gene_Family\tPanaroo_Presence\tPanaroo_Absence'
        header_line += '\tPanP_BlastNoHit\tPanP_BlastHitConcordant\tPanP_BlastHitDiscordant\tPanP_BlastHitMissing\tPanP_BlastOffTarget'
        header_line += '\tPanA_BlastNoHit\tPanA_BlastHitConcordant\tPanA_BlastHitDiscordant\tPanA_BlastHitMissing\tPanA_BlastOffTarget'
        header_line += '\tAvg_Expected_Neighbors\tAvg_Blast_Neighbors\tOverlap_GeneFamilies\tOverlap_GeneFamilies_NoCoOccurance\tDiscordant_Percent_Overlap\n'
        _ = fh.write(header_line)
        for gene_fam, blastdf_gf in blastdf.groupby('Query_SeqID'):
            d = {}
            # count the number of presences and absences in panaroo for this gene family
            gene_fam_row_index = pandf[pandf['Gene'] == gene_fam].index[0]
            pan_data_series = pandf.loc[gene_fam_row_index]
            # remove the Gene column
            pan_data_series = pan_data_series.drop('Gene')
            # get the index names that have non-empty entries - note that panaroo_presence will include refound and pseudo genes
            d['panaroo_presence_isolates'] = set(pan_data_series[pan_data_series != ''].index)
            d['panaroo_absence_isolates'] = set(pan_data_series[pan_data_series == ''].index)
            # iterate though the blastdf_gf and assign each isolate to a category based on the blast hit type
            # nohit means no blast hit at all, concordant means the hit overlaps the expected gene family, discordant means the hit overlaps a different gene family, missing means a hit but no overlap at all
            d['blast_NoHit_isolates'] = set()
            d['blast_HitConcordant_isolates'] = set()
            d['blast_HitDiscordant_isolates'] = set()
            d['blast_HitMissing_isolates'] = set()
            d['blast_OffTarget_isolates'] = set()
            for index, row in blastdf_gf.iterrows():
                isolate_name = row['Subject_SeqID']
                if row['Good_Blast_Hit'] == False:
                    d['blast_NoHit_isolates'].add(isolate_name)
                elif row['Max_Neighbors_in_BLAST'] == 0 and row['Expected_Neighbors'] != 0:
                    d['blast_OffTarget_isolates'].add(isolate_name)
                else:
                    if row['Overlap_Found'] == False:
                        d['blast_HitMissing_isolates'].add(isolate_name)
                    else:
                        if row['Overlap_GeneFamily'] == gene_fam:
                            d['blast_HitConcordant_isolates'].add(isolate_name)
                        else:
                            d['blast_HitDiscordant_isolates'].add(isolate_name)
            if d['panaroo_presence_isolates'] | d['panaroo_absence_isolates'] != d['blast_NoHit_isolates'] | d['blast_HitConcordant_isolates'] | d['blast_HitDiscordant_isolates'] | d['blast_HitMissing_isolates'] | d['blast_OffTarget_isolates']:
                print(f'Error: mismatch in isolate sets for gene family {gene_fam}')
                quit(1)
            out_line = f'{gene_fam}\t{len(d["panaroo_presence_isolates"])}\t{len(d["panaroo_absence_isolates"])}'
            # add the number of NoHit, HitConcordant, etc. isolates that are present in panaroo
            for ikey in ['panaroo_presence_isolates', 'panaroo_absence_isolates']:
                for jkey in ['blast_NoHit_isolates', 'blast_HitConcordant_isolates', 'blast_HitDiscordant_isolates', 'blast_HitMissing_isolates', 'blast_OffTarget_isolates']:
                    count = len(d[ikey] & d[jkey])
                    out_line += f'\t{count}'
            # add average expected and found neighbors
            out_line += f'\t{blastdf_gf["Expected_Neighbors"].mean()}\t{blastdf_gf["Max_Neighbors_in_BLAST"].mean()}'
            # add a string of each unique discordant overlap gene family
            dis_gf = [x for x in blastdf_gf['Overlap_GeneFamily'].unique() if x != gene_fam and x != 'None']
            if dis_gf == []:
                dis_gf = ['NA']
                dis_gf_co = ['NA']
                dis_gf_overlap = -1
            else:
                dis_gf_co = check_discordant_gene_families(gene_fam, dis_gf, pandf)
                # take the mean of the Percent_Overlap column, only counting rows where discordant gene families were found
                dis_gf_df = blastdf_gf[blastdf_gf['Overlap_GeneFamily'].isin(dis_gf)]
                dis_gf_overlap = dis_gf_df['Percent_Overlap'].mean()
            out_line += f'\t{",".join(dis_gf)}\t{",".join(dis_gf_co)}\t{dis_gf_overlap}\n'
            _ = fh.write(out_line)

def check_discordant_gene_families(gene_fam, dis_gf_list, pandf):
    # take a target gene family and a list of discordant gene families
    # for each, check if the target and the dis gf ever co-occur in the same isolate in pandf
    # return a list of the discordant gfs that never co-occur
    out_list = []
    for dis_gf in dis_gf_list:
        if dis_gf == 'None' or dis_gf == '':
            print(f'Error, invalid discordant gene family {dis_gf} for gene family {gene_fam}')
            quit(1)
        target_row = pandf[pandf['Gene'] == gene_fam].iloc[0]
        target_isolates = set(target_row[target_row != ''].index.tolist())
        target_isolates.remove('Gene')
        dis_row = pandf[pandf['Gene'] == dis_gf].iloc[0]
        dis_isolates = set(dis_row[dis_row != ''].index.tolist())
        dis_isolates.remove('Gene')
        if len(target_isolates & dis_isolates) == 0:
            out_list.append(dis_gf)
    if out_list == []:
        out_list = ['NA']
    return out_list

def compare_panaroo_blast_v2(panaroo_file, blast_file, filterfile, output_file):
    # read in the panaroo presence/absence matrix
    pandf = read_panaroo_data(panaroo_file, filterfile)
    # read in the blast results
    blastdf = read_blast_data_v3(blast_file)
    # for each gene family, write out the number of isolates with each type of blast outcome
    with open(output_file, 'w') as out_f:
        header_line_done = False
        for gene_fam in blastdf:
            if gene_fam not in pandf['Gene'].values:
                print(f'{gene_fam} not found in panaroo data')
                quit(1)
            out_data = {}
            gene_fam_row_index = pandf[pandf['Gene'] == gene_fam].index[0]
            pan_data_series = pandf.loc[gene_fam_row_index]
            # remove the Gene column
            pan_data_series = pan_data_series.drop('Gene')
            # get the index names that have non-empty entries - note that panaroo_presence will include refound and pseudo genes
            panaroo_presence = set(pan_data_series[pan_data_series != ''].index)
            panaroo_absence = set(pan_data_series[pan_data_series == ''].index)
            out_data['p_pan'] = len(panaroo_presence)
            out_data['a_pan'] = len(panaroo_absence)
            for key in blastdf[gene_fam]:
                key_presence = len(panaroo_presence & blastdf[gene_fam][key])
                out_data[f'p_{key}'] = key_presence
                key_absence = len(panaroo_absence & blastdf[gene_fam][key])
                out_data[f'a_{key}'] = key_absence
            out_line = ''
            for key in out_data:
                out_line += f'\t{out_data[key]}'
            if not header_line_done:
                header_line = 'Gene_Family'
                for key in out_data:
                    header_line += f'\t{key}'
                out_f.write(f'{header_line}\n')
                header_line_done = True
            out_f.write(f'{gene_fam}{out_line}\n')


def read_blast_data_v3(blast_file):
    # make one dict that describes how many neighbors each blast hit has
    # this is mainly for troubleshooting
    blastdict = {}
    with open(blast_file, 'r') as fhin:
        next(fhin)
        for line in fhin:
            query_seqid, subject_seqid, expected_neighbors, found_neighbors, good_blast_hit, overlap_found, overlap_name, overlap_gf, overlap_category, overlap_percent = line.strip().split('\t')
            if query_seqid not in blastdict:
                blastdict[query_seqid] = {'nohit': set(), '2e2n': set(), '2e1n': set(), '2e0n': set(), '1e1n': set(), '1e0n': set(), '0e0n': set()}
            key = ''
            if good_blast_hit == 'False':
                key = 'nohit'
            else:
                key = f'{expected_neighbors}e{found_neighbors}n'
            blastdict[query_seqid][key].add(subject_seqid)
    return blastdict


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

def check_boolean_arg(input_str):
    if input_str.lower() in ['true', 't', '1']:
        return True
    elif input_str.lower() in ['false', 'f', '0']:
        return False
    else:
        print('Error: Boolean required for --neighbor_mode')
        quit(1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--blast','-b',type=str,
        help='''Provide a BLAST summary file generated by the read_blast script.''',
        required=True
        )
    parser.add_argument(
        '--panaroo','-p',type=str,
        help='''Provide Panaroo presence/absence matrix file in the csv format. This should be the same file the gene families in the 
        input list are drawn from, and the same file used with extract_gene_from_isolate. Gene names should be present in this file.''',
        required=True
        )
    parser.add_argument(
        '--filterfile','-ff',type=str,
        help='''Provide a path to a one-column csv consisting only of isolate names. Only isolates in this file will be used.''',
        default=None
        )
    parser.add_argument(
        '--output_overlap','-oo',type=str,
        help='''Provide a file name for the output file for the overlap table.''',
        required=True
        )
    parser.add_argument(
        '--output_neighbors','-on',type=str,
        help='''Provide a file name for the output file for the neighbor table.''',
        required=True
        )
    args = parser.parse_args()
    compare_panaroo_blast_v3(args.panaroo, args.blast, args.filterfile, args.output_overlap)
    compare_panaroo_blast_v2(args.panaroo, args.blast, args.filterfile, args.output_neighbors)


if __name__ == '__main__':
    main()