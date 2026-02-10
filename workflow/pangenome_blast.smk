
configfile: "config/config.yaml"

import pandas as pd
import os

PREFIX = config["prefix"]
CLADE = config["clade"]

rule all:
    input:
        comparison_output_high = expand("results/{PREFIX}/output2/panaroo_blast_comparison_{PREFIX}_{CLADE}_high_confidence.tsv", PREFIX=PREFIX, CLADE=CLADE),
        #comparison_output_med = expand("results/{PREFIX}/output2/panaroo_blast_comparison_MediumConfidence_{PREFIX}_{CLADE}.tsv", PREFIX=PREFIX, CLADE=CLADE),
        #blast_summary_output = expand("results/{PREFIX}/output1/blast_summary_output_{PREFIX}_{CLADE}.fasta", PREFIX=PREFIX, CLADE=CLADE),
        #pangenome_neighbors_complete = expand("results/{PREFIX}/neighbor_db/pangenome_neighbors_complete_{PREFIX}_{CLADE}.txt", PREFIX=PREFIX, CLADE=CLADE),
        #comparison_output = expand("results/{PREFIX}/output2/panaroo_blast_comparison_{PREFIX}_{CLADE}.fasta", PREFIX=PREFIX, CLADE=CLADE),


rule extract_accessory_gene_sequences:
    input:
        accessory_gene_list = config["accessory_gene_list"],
    output:
        accessory_gene_sequences = "results/{PREFIX}/output1/extracted_accessory_sequences_{PREFIX}_{CLADE}.fasta",
    params:
        assemblies_dir = config["assemblies_dir"],
        gff_dir = config["gff_dir"],
        matrix = config["pangenome_matrix"],
        filterfile = config.get("filterfile", None),
    resources:
        mem_mb = 4000,
        runtime = 60,
    threads: 1
    conda:
        'biopython'
    shell:
        """
        python3.12 scripts/extract_gene_from_isolate.py --input {input.accessory_gene_list} --output {output.accessory_gene_sequences} \
        --assemblies {params.assemblies_dir} --gff {params.gff_dir} --pangenome {params.matrix} --filterfile {params.filterfile}
        """


rule download_caurisblast:
    output:
        blast_script = "caurisblast/blast.py",
    params:
        caurisblast_repo = config["caurisblast_repo"],
    resources:
        mem_mb = 2000,
        runtime = 10,
    threads: 1
    conda:
        'biopython'
    shell:
        "mkdir -p caurisblast && git clone {params.caurisblast_repo} caurisblast"


rule caurisblast:
    input:
        accessory_gene_sequences = "results/{PREFIX}/output1/extracted_accessory_sequences_{PREFIX}_{CLADE}.fasta",
        blast_script = "caurisblast/blast.py",
    output:
        blast_raw_output = "results/{PREFIX}/output1/blast_raw_output_{PREFIX}_{CLADE}.tsv",
        #blast_raw_output = "caurisblast/results/{PREFIX}/{PREFIX}_nucl_blastn_1e-5_blast_results.csv",
    params:
        blast_subject = config["blast_subject"],
        eval_threshold = config["eval_threshold"],
        wordsize = config["wordsize"],
        filterfile = config.get("filterfile", None),
        batch_name = f"{PREFIX}_{CLADE}",
    resources:
        mem_mb = 6000,
        runtime = 120,
    conda:
        'biopython'
    threads: 4
    shell:
        """
        # if the subject is a fasta file, use it with the -d flag 
        # otherwise, assume blastdb generation is needed and exclude the -d flag
        # caurisblast/results/{params.batch_name}/{params.batch_name}_nucl_blastn_{params.eval_threshold}_blast_results.csv
        if [[ "{params.blast_subject}" == *.fasta || "{params.blast_subject}" == *.fa ]]; then
            python3.12 caurisblast/blast.py -q {input.accessory_gene_sequences} -s {params.blast_subject} \
            -th {threads} -e {params.eval_threshold} -w {params.wordsize} -k blastn -t nucl -v \
            -n {params.batch_name} -d
        else
            python3.12 caurisblast/blast.py -q {input.accessory_gene_sequences} -s {params.blast_subject} \
            -th {threads} -e {params.eval_threshold} -w {params.wordsize} -k blastn -t nucl -v \
            -n {params.batch_name} --filterfile {params.filterfile}
        fi
        mv caurisblast/results/{params.batch_name}/{params.batch_name}_nucl_blastn_{params.eval_threshold}_blast_results.csv {output.blast_raw_output}
        """

rule find_pangenome_neighbors:
    input:
        blast_raw_output = "results/{PREFIX}/output1/blast_raw_output_{PREFIX}_{CLADE}.tsv",
    output:
        pangenome_neighbors_complete = "results/{PREFIX}/neighbor_db/{CLADE}/pangenome_neighbors_complete_{PREFIX}_{CLADE}.txt",
    params:
        out_dir = "results/{PREFIX}/neighbor_db/{CLADE}/",
        accessory_gene_list = config["accessory_gene_list"],
        isolate_list = config["filterfile"],
        pan_matrix = config["pangenome_matrix"],
        pan_graph = config["pangenome_graph"],
        gff_dir = config["gff_dir"],
    resources:
        mem_mb = 40000,
        runtime = 600,
    threads: 1
    conda:
        'biopython'
    shell:
        """
        python3.12 scripts/find_pangenome_neighbors.py --accessory_list {params.accessory_gene_list} --isolate_list {params.isolate_list} \
        --pan_matrix {params.pan_matrix} --pan_graph {params.pan_graph} --gff_dir {params.gff_dir} \
        --out_dir {params.out_dir}
        exitcode=$?
        if [ $exitcode == 0 ]; then
            touch {output.pangenome_neighbors_complete}
        else
            echo "Error in finding pangenome neighbors"
            exit $exitcode
        fi
        """

rule summarize_neighbor_db:
    input:
        pangenome_neighbors_complete = "results/{PREFIX}/neighbor_db/{CLADE}/pangenome_neighbors_complete_{PREFIX}_{CLADE}.txt",
    output:
        neighbor_summary = "results/{PREFIX}/output2/neighborDB_summary_{PREFIX}_{CLADE}.tsv",
    params:
        neighbor_db = "results/{PREFIX}/neighbor_db/{CLADE}/",
    resources:
        mem_mb = 15000,
        runtime = 600,
    threads: 1
    conda:
        'biopython'
    shell:
        """
        python3.12 scripts/summarize_neighbor_db.py --neighbor_db {params.neighbor_db} --output {output.neighbor_summary}
        """

rule summarize_blast_results:
    input:
        blast_raw_output = "results/{PREFIX}/output1/blast_raw_output_{PREFIX}_{CLADE}.tsv",
        neighbor_summary = "results/{PREFIX}/output2/neighborDB_summary_{PREFIX}_{CLADE}.tsv",
        pangenome_neighbors_complete = "results/{PREFIX}/neighbor_db/{CLADE}/pangenome_neighbors_complete_{PREFIX}_{CLADE}.txt",
    output:
        blast_summary_output = "results/{PREFIX}/output1/blast_summary_output_{PREFIX}_{CLADE}.tsv",
    params:
        neighbor_db = "results/{PREFIX}/neighbor_db/{CLADE}/",
        minimum_identity = config.get("minimum_identity", 0.9),
        minimum_evalue = config.get("minimum_evalue", 1e-5),
        minimum_coverage = config.get("minimum_coverage", 0.9),
        max_distance = config.get("max_distance", 50000),
    resources:
        mem_mb = 15000,
        runtime = 600,
    threads: 1
    conda:
        'biopython'
    shell:
        """
        python3.12 scripts/read_blast_neighbors_v2.py --input {input.blast_raw_output} --output {output.blast_summary_output} \
        --neighbor_db {params.neighbor_db} --minimum_identity {params.minimum_identity} --minimum_evalue {params.minimum_evalue} \
        --minimum_coverage {params.minimum_coverage} --maximum_distance {params.max_distance}
        """


rule compare_panaroo_blast:
    input:
        blast_summary_output = "results/{PREFIX}/output1/blast_summary_output_{PREFIX}_{CLADE}.tsv",
    output:
        comparison_output_high = "results/{PREFIX}/output2/panaroo_blast_comparison_{PREFIX}_{CLADE}_high_confidence.tsv",
        comparison_output_medium = "results/{PREFIX}/output2/panaroo_blast_comparison_{PREFIX}_{CLADE}_medium_confidence.tsv",
        comparison_output_low = "results/{PREFIX}/output2/panaroo_blast_comparison_{PREFIX}_{CLADE}_low_confidence.tsv",
    params:
        matrix = config["pangenome_matrix"],
        filterfile = config.get("filterfile", None),
        neighbor_mode = config.get("neighbor_mode", False),
        confidence_level = 'all',
        output_filename = "results/{PREFIX}/output2/panaroo_blast_comparison_{PREFIX}_{CLADE}.tsv",
    resources:
        mem_mb = 8000,
        runtime = 600,
    threads: 1
    conda:
        'biopython'
    shell:
        """
        python3.12 scripts/compare_panaroo_blast_v2.py --blast {input.blast_summary_output} --panaroo {params.matrix} \
        --output {params.output_filename} --filterfile {params.filterfile} --neighbor_mode {params.neighbor_mode} \
        --confidence {params.confidence_level}
        """


