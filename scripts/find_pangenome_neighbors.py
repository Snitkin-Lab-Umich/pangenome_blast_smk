
import argparse
import pandas as pd
from itertools import filterfalse
import networkx as nx
import gffutils as gff


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

def gene_name_checker(gene_name):
    # determine if the input string is one that is present in the original gff
    # our assemblies append '.cds' to the end of each gene
    # the b8441 reference appends '.1' instead
    # if a valid gene is found, return the gene name, otherwise return False
    if gene_name == '':
        return False
    if ';' in gene_name:
        gene_list = gene_name.split(';')
    else:
        gene_list = [gene_name]
    # remove the _pseudo suffix - pseudogenes still have entries in the gff and will function as neighboring genes if present
    gene_list = [x.replace('_pseudo', '') for x in gene_list]
    present_gene_list = [x for x in gene_list if x.endswith(('.cds','.1'))]
    if present_gene_list == []:
        return False
    else:
        # if there are multiple valid genes in this gene family (i.e. paralogs), simply take the first one
        # (this could give misleading results if paralogs are on different chromosomes)
        return present_gene_list[0]

def gene_name_checker_paralogskip(gene_name):
    # determine if the input string is one that is present in the original gff
    # our assemblies append '.cds' to the end of each gene
    # the b8441 reference appends '.1' instead
    # if a valid gene is found, return the gene name, otherwise return False
    # skip any gene families with multiple paralogs, treating them as absences
    if gene_name == '' or ';' in gene_name:
        return False
    # remove the _pseudo suffix - pseudogenes still have entries in the gff and will function as neighboring genes if present
    gene_name = gene_name.replace('_pseudo', '')
    if gene_name.endswith(('.cds','.1')):
        return gene_name
    else:
        return False

def neighbor_search(gene_family, isolate, neighbor_dict, gff_db, pan_matrix):
    # get the neighboring gene families for this query gene family
    neighbor_list = neighbor_dict.get(gene_family, [])
    #print(f'initial neighbor list: {neighbor_list}')
    present_neighbor_dict = {}
    for neighbor in neighbor_list:
        # print(f'finding present neighbor for gene family {gene_family} neighbor {neighbor} in isolate {isolate}')
        # check if these neighbors are actually present in this assembly, based on the pangenome matrix
        # if they are not, attempt to locate the nearest present neighbor by traversing the pangenome graph
        # start with a list of genes to check for presence
        target_genefam_list = [neighbor]
        target_genefam = None
        # make sure we do not traverse backwards towards the original gene family
        genes_seen = [gene_family]
        target_genefam_GeneName = False
        while not target_genefam_GeneName:
            # if there are no more neighbors to check, or if we have searched more than 20 genes, give up and mark as not_found
            if target_genefam_list == [] or len(genes_seen) > 20:
                present_neighbor_dict[neighbor] = 'not_found'
                target_genefam = None
                break
            # take the first gene family from the list
            target_genefam = target_genefam_list.pop(0)
            # add to the list of genes already seen
            genes_seen.append(target_genefam)
            # check if this gene family is present
            target_genefam_GeneName_string = pan_matrix.loc[pan_matrix['Gene'] == target_genefam, isolate].iloc[0]
            target_genefam_GeneName = gene_name_checker_paralogskip(target_genefam_GeneName_string)
            # find the next neighbors of this gene family (which can be none or multiple)
            next_neighbors = neighbor_dict.get(target_genefam, [])
            # exclude any previously seen genes
            next_neighbors = [x for x in next_neighbors if x not in genes_seen]
            for n in next_neighbors:
                if n not in target_genefam_list:
                    target_genefam_list.append(n)
        # if a present neighbor was successfully found, add it to the present_neighbor_dict
        # note that it is possible for multiple neighbors to resolve into the same target_genefam through this method (if the pangenome graph makes a loop)
        # if this happens, repeatedly append _REP to the gene name until it is unique
        if target_genefam is not None:
            while target_genefam in present_neighbor_dict.keys():
                target_genefam += '_REP'
            present_neighbor_dict[target_genefam] = target_genefam_GeneName
    neighbor_data = resolve_graph_neighbors(neighbor_list, present_neighbor_dict, gff_db, isolate)
    return neighbor_data


def resolve_graph_neighbors(graph_neighbor_list, found_neighbor_dict, gff_db, isolate):
    # take a list of neighboring gene families from the pangenome graph
    # compare this to a dictionary of neighboring gene names
    # return a list of neighboring gene data (gene family, gene name, scaffold, start, end) for exactly two genes
    # distinguish between cases where there are no neighbors in the graph list vs cases where a present neighbor was unable to be found
    found_neighbor_data = []
    for neighbor_GeneFam, neighbor_GeneName in found_neighbor_dict.items():
        if neighbor_GeneName == 'not_found':
            found_neighbor_data.append((neighbor_GeneFam,'not_found','absent','absent','absent'))
        else:
            gene_feature = gff_db[neighbor_GeneName]
            scaffold = gene_feature.seqid
            start = gene_feature.start
            end = gene_feature.end
            found_neighbor_data.append((neighbor_GeneFam, neighbor_GeneName, scaffold, start, end))
    # if there are at least two found neighbors, found_neighbor_data can be returned without adding anything
    # if there are two neighbors that share a scaffold, try to return those two
    if len(graph_neighbor_list) != len(found_neighbor_data):
        print(f'Error: mismatch between pangenome graph neighbors and located neighbors for isolate {isolate}.')
    if len(found_neighbor_data) >= 2:
        # first, attempt to return any two neighbors that share a scaffold
        for i in range(len(found_neighbor_data)):
            for j in range(i+1, len(found_neighbor_data)):
                if found_neighbor_data[i][2] == found_neighbor_data[j][2] and found_neighbor_data[i][1] != 'not_found':
                    return [found_neighbor_data[i], found_neighbor_data[j]]
        # next, move any neighbors that were not found to the end of the list
        found_neighbor_data.sort(key=lambda x: x[1] == 'not_found')
        # return the first two neighbors 
        return found_neighbor_data[:2]
    # if there are still less than two neighbors, it means the pangenome graph had less than two neighbors
    while len(found_neighbor_data) < 2:
        found_neighbor_data.append(('no_neighbor','no_neighbor','absent','absent','absent'))
    return found_neighbor_data[:2]


def create_neighbor_table(accessory_gene_list, isolate_list, pan_matrix, pan_graph, gff_dir, output_dir, spliced_gff):
    # for each isolate in the isolate_list, iterate through the full list of accessory gene families
    # for each gene family, identify the position of the gene in this isolate (if present)
    # then, identify the positions of the neighboring gene families as well (if present)
    # finally, output a file containing the positions of the target gene, both of its neighbors, and if both neighbors are on the same contig/scaffold
    ## first, get the pangenome matrix and pangenome graph
    neighbor_dict = read_pangenome_graph(pan_graph)
    pan_matrix = read_panaroo_data(pan_matrix)
    for isolate in isolate_list:
        output_file = f'{output_dir}/{isolate}_accessory_neighborfile.tsv'
        with open(output_file, 'w') as fhout:
            fhout.write('GeneFamily\tGeneName\tScaffold\tStart\tEnd\tNeighbor1_GeneFamily\tNeighbor1_GeneName\tNeighbor1_Scaffold\tNeighbor1_Start\tNeighbor1_End\tNeighbor2_GeneFamily\tNeighbor2_GeneName\tNeighbor2_Scaffold\tNeighbor2_Start\tNeighbor2_End\tNeighbor_Shared_Scaffold\tPresent_Shared_Scaffold\n')
            # get the gffdb for this isolate
            gff_file = f'{gff_dir}/{isolate}.gff'
            if not spliced_gff:
                gff_file = f'{gff_dir}/{isolate}.gff3'
            gff_db = gff.create_db(gff_file,dbfn=":memory:",force=True,keep_order=False,merge_strategy="create_unique",sort_attribute_values=True,from_string=False)
            for genefam in accessory_gene_list:
                # determine if this gene family is present in this isolate
                genefam_GeneName_string = pan_matrix.loc[pan_matrix['Gene'] == genefam, isolate].iloc[0]
                genefam_GeneName = gene_name_checker(genefam_GeneName_string)
                if not genefam_GeneName:
                    genefam_data = (genefam,'absent','absent','absent','absent')
                else:
                    gene_feature = gff_db[genefam_GeneName]
                    scaffold = gene_feature.seqid
                    start = gene_feature.start
                    end = gene_feature.end
                    genefam_data = (genefam, genefam_GeneName, scaffold, start, end)
                # regardless of presence/absence of the target gene, determine the position of the neighbors
                neighbor_data = neighbor_search(genefam, isolate, neighbor_dict, gff_db, pan_matrix)
                # neighbor data should now always contain exactly two entries
                neighbor_shared_scaffold = 'FALSE'
                if neighbor_data[0][2] == neighbor_data[1][2] and neighbor_data[0][2] != 'absent':
                    neighbor_shared_scaffold = 'TRUE'
                present_shared_scaffold = 'FALSE'
                if genefam_data[2] == neighbor_data[0][2] and genefam_data[2] == neighbor_data[1][2] and genefam_data[2] != 'absent':
                    present_shared_scaffold = 'TRUE'
                outdata = genefam_data + neighbor_data[0] + neighbor_data[1] + (neighbor_shared_scaffold,present_shared_scaffold)
                outdata = [str(x) for x in outdata]
                outstring = '\t'.join(outdata) + '\n'
                fhout.write(outstring)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--accessory_list','-a',type=str,
        help='''Provide a one-column csv consisting of the accessory gene families to search for. This should have a header "gene_family".''',
        default=None
        )
    parser.add_argument(
        '--isolate_list','-i',type=str,
        help='''Provide a one-column csv consiting of the isolates to search within. This should have no header.''',
        default=None
        )
    parser.add_argument(
        '--pan_matrix','-pm',type=str,
        help='''Provide a path to a pangenome presence/absence matrix file. Gene names should be present.''',
        default=None
        )
    parser.add_argument(
        '--pan_graph','-pg',type=str,
        help='''Provide a path to a pangenome graph file in gml format.''',
        default=None
        )
    parser.add_argument(
        '--gff_dir','-gff',type=str,
        help='''Provide a path to a directory of GFF files used in the generation of the pangenome data.''',
        default=None
        )
    parser.add_argument(
        '--out_dir','-o',type=str,
        help='''Provide a path to a directory for output files.''',
        default=None
        )
    parser.add_argument(
        '--spliced_gff','-sg',type=str,choices=['True','False'],
        help='''Provide a boolean value to indicate whether to use spliced GFF files.''',
        default='True'
        )
    args = parser.parse_args()
    if args.spliced_gff == 'True':
        args.spliced_gff = True
    else:
        args.spliced_gff = False
    # read in the accessory gene list
    accessory_df = pd.read_csv(args.accessory_list)
    accessory_gene_list = accessory_df['gene_family'].tolist()
    # read in the isolate list
    isolate_df = pd.read_csv(args.isolate_list,header=None)
    isolate_list = isolate_df[0].tolist()
    # create the neighbor table
    create_neighbor_table(accessory_gene_list, isolate_list, args.pan_matrix, args.pan_graph, args.gff_dir, args.out_dir, args.spliced_gff)

if __name__ == '__main__':
    main()
