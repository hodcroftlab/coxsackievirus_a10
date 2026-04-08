###############
# Snakemake execution templates:

# To run a default protein xy run:
# snakemake  auspice/coxsackievirus_A10_vp1.json --cores 9

# To run a default whole genome run (>6400bp):
# snakemake auspice/coxsackievirus_A10_genome.json --cores 9

import os
from datetime import date

# Load config file
if not config:
    configfile: "config/config.yaml"

# Load environment variables
# Try to load .env, but don't fail if it doesn't exist (for Actions)
try:
    from dotenv import load_dotenv
    load_dotenv(".env")
except:
    pass

REMOTE_GROUP = os.getenv("REMOTE_GROUP")
UPLOAD_DATE = date.today().isoformat()

FETCH_SEQUENCES=True

###############
#ensure vp1 name similar to that found in the reference_sequence.gb CDS
wildcard_constraints:
    seg="vp1|whole_genome"  # Define segments to analyze, e.g. vp1, whole-genome. This wildcard will be used in the rules "{seg}" to define the path or protein to use
   
# Define segments to analyze
segments = ['vp1', 'whole-genome'] # This is only for the expand in rule all

# Rule to handle configuration files and data file paths
rule files:
    input:
        sequence_length =   "{seg}",
        colors =            "config/colors.tsv",
        dropped_strains =   "config/dropped_strains.txt",
        regions=            "config/geo_regions.tsv",
        lat_longs =         "config/lat_longs.tsv",
        reference =         "{seg}/config/reference_sequence.gb", 
        gff_reference =     "{seg}/config/annotation.gff3",
        auspice_config =    "{seg}/config/auspice_config.json",
        clades =            "{seg}/config/clades_genome.tsv",
        include =           "config/include.txt",
        SEQUENCES =         "data/sequences.fasta",
        METADATA =          "data/metadata.tsv",
        meta_collab =       "data/meta_collab.tsv", # collaborator metadata, empty for now
        meta_publications = "data/meta_publications.tsv", # Published metadata
        last_updated_file = "data/date_last_updated.txt",
        local_accn_file =   "data/local_accn.txt",

files = rules.files.input

# Expand augur JSON paths
rule all:
    input:
        augur_jsons = expand("auspice/coxsackievirus_A10_{segs}.json", segs=segments),
        meta = files.METADATA,
        seq = files.SEQUENCES

##############################
# Download from NBCI Virus with ingest snakefile
###############################
if FETCH_SEQUENCES == True:
    rule fetch:
        input:
            dir = "ingest"
        output:
            sequences=files.SEQUENCES,
            metadata=files.METADATA
        threads: workflow.cores
        shell:
            """
            cd {input.dir} 
            snakemake --cores {threads} all
            cd ../
            """

##############################
# Optional: Fetch metadata from genbank
###############################

# This rule is very slow. Only give accessions as input where you are certain that they have GenBank metadata.
# rule fetch_metadata:
#     message:
#         """
#         Retrieving GenBank metadata for the specified accessions.
#         """
#     input:
#         accessions="data/metadata/afm.txt",
#         config="config/config.yaml" # include symptom list and isolation source mapping
#     output:
#         metadata="data/metadata/afm.tsv",
#     params:
#         virus="Coxsackievirus A10",
#         genbank_metadata="data/genbank_metadata.tsv",
#         columns = "accession strain published_clade country location date age gender diagnosis doi"

#     log:
#         "logs/fetch_metadata.log"
#     shell:
#         """
#         python scripts/fetch_genbank_metadata.py \
#             --virus "{params.virus}" \
#             --accession_file {input.accessions} \
#             --output {output.metadata} \
#             --genbank {params.genbank_metadata} \
#             --config {input.config} \
#             --columns {params.columns} \
#             2> {log}
#         """

##############################
# AUGUR CURATE AND MERGE
# Change the format of the dates in the metadata
# Attention: ```augur curate``` only accepts iso 8 formats; please make sure that you save e.g. Excel files in the correct format
# Merge with other metadata files you might have
###############################

rule curate:
    message:
        """
        Cleaning up metadata with augur curate
        """
    input:
        metadata = files.METADATA,  # Path to input metadata file
        meta_publications = files.meta_publications,
        genbank_metadata="data/genbank_metadata.tsv"
    params:
        strain_id_field=config["id_field"],
        date_fields=config["curate"]["date_fields"],
        expected_date_formats=config["curate"]["expected_date_formats"],
        temp=temp("temp/merged_metadata.tsv")
    output:
        metadata="data/merged_meta.tsv"
    shell:
        """               
        # Merge curated metadata
        augur merge --metadata metadata={input.metadata} meta_publications={input.meta_publications} genbank={input.genbank_metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata {params.temp}

        # Normalize strings for publication metadata
        augur curate normalize-strings \
            --id-column {params.strain_id_field} \
            --metadata {params.temp} \
        | augur curate format-dates \
            --date-fields {params.date_fields} \
            --no-mask-failure \
            --expected-date-formats {params.expected_date_formats} \
            --id-column {params.strain_id_field} \
            --output-metadata {output.metadata}
        """

##############################
# Add additional sequences
# if you have sequences that are not on NCBI Virus
###############################

rule update_sequences:
    input:
        sequences = files.SEQUENCES,
        metadata= rules.curate.output.metadata,
    output:
        sequences = "data/all_sequences.fasta",
        metadata = "data/all_metadata.tsv"
    params:
        strain_id_field=config["id_field"],
        file_ending = "data/*.fas*",
        temp = temp("temp/sequences.fasta"),
        date_last_updated = files.last_updated_file,
        local_accn = files.local_accn_file,
    shell:
        """
        touch {params.temp} && rm {params.temp}
        cat {params.file_ending} > {params.temp}
        python scripts/update_sequences.py --in_seq {params.temp} --out_seq {output.sequences} --dates {params.date_last_updated} \
        --local_accession {params.local_accn} --meta {input.metadata} --ingest_seqs {input.sequences}

        awk '/^>/{{if (seen[$1]++ == 0) print; next}} !/^>/{{print}}' {output.sequences} > {params.temp} && mv {params.temp} {output.sequences}

        augur merge --metadata metadata={input.metadata} dates={params.date_last_updated} \
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata {output.metadata}
        """

##############################
# BLAST
# blast fasta files for your specific proteins
# cut out your protein from fasta sequences
###############################

rule extract:
    input: 
        genbank_file = files.reference
    output: 
        extracted_fasta = "{seg}/results/extracted.fasta"    
    params:
        product_name = "{seg}"
    shell:
        """
        python scripts/extract_gene_from_whole_genome.py \
        --genbank_file {input.genbank_file} \
        --output_fasta {output.extracted_fasta} \
        --product_name {params.product_name}

        """

rule blast:
    input: 
        blast_db_file = rules.extract.output.extracted_fasta,  
        seqs_to_blast = rules.update_sequences.output.sequences
    output:
        blast_out = "temp/{seg}/blast_out.csv"
    params:
        blast_db = "temp/{seg}/blast_database"
    shell:
        """
        sed -i 's/-//g' {input.seqs_to_blast}
        makeblastdb -in {input.blast_db_file} -out {params.blast_db} -dbtype nucl
        blastn -task blastn -query {input.seqs_to_blast} -db {params.blast_db}\
        -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart \
        send evalue bitscore qcovs' -out {output.blast_out} -evalue 0.0005
        """

rule blast_sort: #TODO: change the parameters in blast_sort.py (replace lengths with your specific protein)
    input:
        blast_result = rules.blast.output.blast_out, # output blast (for your protein)
        seqs_to_blast = rules.update_sequences.output.sequences
    output:
        sequences = "{seg}/results/sequences.fasta"
    params:
        range = "{seg}",  # Determines which protein (or whole genome) is processed
        min_length = lambda wildcards: {"vp1": 600, "whole_genome": 6400}[wildcards.seg],  # Min length ## DONE replace whole_genome with genome, min length = 6000
        max_length = lambda wildcards: {"vp1": 900, "whole_genome": 8000}[wildcards.seg]    
    shell:
        """
        python scripts/blast_sort.py --blast {input.blast_result} \
            --seqs {input.seqs_to_blast} \
            --out_seqs {output.sequences} \
            --range {params.range} \
            --min_length {params.min_length} \
            --max_length {params.max_length}
        """

##############################
# Indexing sequences and filter them.
###############################


rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering
        """
    input:
        sequences = rules.blast_sort.output.sequences
    output:
        sequence_index = "{seg}/results/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - from {params.min_date} onwards
          - excluding strains in {input.exclude}
        """
    input:
        sequences = rules.blast_sort.output.sequences,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = rules.update_sequences.output.metadata,
        exclude = files.dropped_strains,
        include = files.include,
    output:
        sequences = "{seg}/results/filtered.fasta",
        reason ="{seg}/results/filter_log.tsv"
    params:
        group_by = "country year",
        sequences_per_group = 200, # add a limit per group
        strain_id_field= config["id_field"],
        min_date = 1980,  # add a reasonable min date
        min_length = lambda wildcards: {"vp1": 600, "whole_genome": 6400}[wildcards.seg], 
        max_length = lambda wildcards: {"vp1": 900, "whole_genome": 8000}[wildcards.seg]  
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --exclude {input.exclude} \
            --include {input.include} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --min-length {params.min_length} --max-length {params.max_length} \
            --output-sequences {output.sequences} \
            --output-log {output.reason}
        """
# --exclude-where ... or other parameters can be added, see `augur filter --h` for more options

##############################
# Reference for alignment added to sub-folders
###############################

rule reference_gb_to_fasta:
    message:
        """
        Converting reference sequence from genbank to fasta format and putting it in the reference folders of your proteins
        """
    input:
        reference = files.reference

    output:
        reference = "{seg}/results/reference_sequence.fasta"
    run:
        from Bio import SeqIO
        SeqIO.convert(input.reference, "genbank", output.reference, "fasta")

rule align: 
    message:
        """
        Aligning sequences to {input.reference} using Nextclade run.
        """
    input:
        gff_reference = files.gff_reference,
        sequences = rules.filter.output.sequences,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        alignment = "{seg}/results/aligned.fasta",
        tsv = "{seg}/results/nextclade.tsv",    
    params:
        penalty_gap_extend = config["align"]["penalty_gap_extend"],
        penalty_gap_open = config["align"]["penalty_gap_open"],
        penalty_gap_open_in_frame = config["align"]["penalty_gap_open_in_frame"],
        penalty_gap_open_out_of_frame = config["align"]["penalty_gap_open_out_of_frame"],
        kmer_length = config["align"]["kmer_length"],
        kmer_distance = config["align"]["kmer_distance"],
        min_match_length = config["align"]["min_match_length"],
        allowed_mismatches = config["align"]["allowed_mismatches"],
        min_length = config["align"]["min_length"]
    threads: workflow.cores
    shell:
        """
        nextclade3 run \
        -j {threads} \
        {input.sequences} \
        --input-ref {input.reference} \
        --input-annotation {input.gff_reference} \
        --penalty-gap-open {params.penalty_gap_open} \
        --penalty-gap-extend {params.penalty_gap_extend} \
        --penalty-gap-open-in-frame {params.penalty_gap_open_in_frame} \
        --penalty-gap-open-out-of-frame {params.penalty_gap_open_out_of_frame} \
        --kmer-length {params.kmer_length} \
        --kmer-distance {params.kmer_distance} \
        --min-match-length {params.min_match_length} \
        --allowed-mismatches {params.allowed_mismatches} \
        --min-length {params.min_length} \
        --include-reference false \
        --output-tsv {output.tsv} \
        --output-translations "{wildcards.seg}/results/translations/cds_{{cds}}.translation.fasta" \
        --output-fasta {output.alignment}
        """


##############################
# Building a tree
###############################

rule tree:
    message:
        """
        Creating a maximum likelihood tree
        """
    input:
        alignment = rules.align.output.alignment

    output:
        tree = "{seg}/results/tree_raw.nwk"

    threads: workflow.cores
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --nthreads {threads}\
            --output {output.tree}
        """

##############################
# Refine to a timeline
###############################

rule refine:
    message:
        """
        Refining tree by rerooting and resolving polytomies
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output.alignment,
        metadata= rules.update_sequences.output.metadata,
    output:
        tree = "{seg}/results/tree.nwk",
        node_data = "{seg}/results/branch_lengths.json"
    params:
        coalescent = "opt",
        rooting = "mid_point",  # or use a specific accession ID
        date_inference = "marginal",
        clock_filter_iqd = lambda wildcards: {"vp1": 4, "whole_genome": 8}[wildcards.seg],  # set to 6 if you want more control over outliers
        strain_id_field = config["id_field"],
        clock_rate = 0.004, # clockor2 (2.7–4.2e-3)
        clock_std_dev = 0.0015

    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root {params.rooting} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --clock-rate {params.clock_rate}\
            --clock-std-dev {params.clock_std_dev} \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

##############################
# Ancestral sequences and amino acids
###############################

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment

    output:
        node_data = "{seg}/results/nt_muts.json"

    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --keep-ambiguous\
            --inference {params.inference}
        """
 
rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "{seg}/results/aa_muts.json"

    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data}
        """

rule traits:
    message: "Inferring ancestral traits for {params.traits!s}"
    input:
        tree = rules.refine.output.tree,
        metadata= rules.update_sequences.output.metadata
    output:
        node_data = "{seg}/results/traits.json"
        
    params:
        traits = "country",
        strain_id_field= config["id_field"]
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-node-data {output.node_data} \
            --columns {params.traits} \
            --confidence
        """

##############################
# Assign clades or subgenotypes based on list provided
###############################
rule clades: 
    message: "Assigning clades according to nucleotide mutations"
    input:
        tree=rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = files.clades # TODO: assign mutations to specific clades
    output:
        clade_data = "{seg}/results/clades.json"

    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
        """

rule add_url_to_metadata:
    input:
        rules.update_sequences.output.metadata
    output:
        "data/final_meta_url.tsv"
    run:
        import pandas as pd

        meta = pd.read_csv(input[0], sep="\t", dtype=str)
        meta['url'] = "https://www.ncbi.nlm.nih.gov/nuccore/" + meta['accession']
        meta.to_csv(output[0], sep="\t", index=False)

#########################
#  EXPORT
#########################
rule export:
    message: "Creating auspice JSONs"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.add_url_to_metadata.output,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        clades = rules.clades.output.clade_data,
        colors = files.colors,
        lat_longs = files.lat_longs,
        auspice_config = files.auspice_config
    params:
        strain_id_field = config["id_field"],
    output:
        auspice_json = "auspice/coxsackievirus_A10_{seg}.json"
        
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} \
                {input.aa_muts} {input.clades} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """

rule rename_whole_genome:
    message: "Rename whole-genome built"
    input: 
        json="auspice/coxsackievirus_A10_whole_genome.json"
    output:
        json="auspice/coxsackievirus_A10_whole-genome.json" # easier view in auspice
    shell:
        """
        mv {input.json} {output.json}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "ingest/data/*.*",
        "*/results/*",
        "auspice/*.json",
        "temp/*",
        "logs/*",
        "benchmark/*",
        files.METADATA,
        files.SEQUENCES,
        "data/curated/*",
        "data/all_sequences.fasta",
        "data/all_metadata.tsv",
        "data/merged_metadata.tsv",
        "logs/*"
    shell:
        "rm -rfv {params}"

rule upload: ## make sure you're logged in to Nextstrain
    message: "Uploading auspice JSONs to Nextstrain"
    input:
        jsons = expand("auspice/coxsackievirus_A16_{segs}.json", segs=segments)
    params:
        remote_group=REMOTE_GROUP,
        date=UPLOAD_DATE,
        USERNAME=os.getenv("NEXTSTRAIN_REMOTE_USERNAME"),

    shell:
        """
        nextstrain login --username {params.USERNAME}
        nextstrain remote upload \
            nextstrain.org/groups/{params.remote_group}/ \
            {input.jsons}
        nextstrain logout
        mkdir -p auspice/{params.date}
        cp {input.jsons} auspice/{params.date}/
        """