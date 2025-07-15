###############
# Snakemake execution templates:

# To run a default protein xy run:
# snakemake  auspice/<your_virus>_protein_xy.json --cores 9

# To run a default whole genome run (>6400bp):
# snakemake auspice/<your_virus>_genome.json --cores 9

# Load config file
if not config:
    configfile: "config/config.yaml"

###############
#ensure protein_xy name similar to that found in the reference_sequence.gb CDS
wildcard_constraints:
    seg="vp1|whole_genome"  # Define segments to analyze, e.g. vp1, whole-genome. This wildcard will be used in the rules "{seg}" to define the path or protein to use
   
# Define segments to analyze
segments = ['vp1', 'whole_genome'] # This is only for the expand in rule all

# Expand augur JSON paths
rule all:
    input:
        sequences="data/sequences.fasta",
        metadata="data/metadata.tsv",
        augur_jsons = expand("auspice/cva10_{segs}.json", segs=segments) ## TODO: replace <your_virus> with actual virus name (Ctrl+H)


# Rule to handle configuration files and data file paths
rule files:
    input:
        sequence_length =       "{seg}",
        colors =                "config/colors.tsv",
        dropped_strains =       "config/dropped_strains.txt",
        regions=                "config/geo_regions.tsv",
        lat_longs =             "config/lat_longs.tsv",
        reference =             "{seg}/config/reference_sequence.gb", ####TODO: provide a reference sequence
        auspice_config =        "{seg}/config/auspice_config.json",
        clades =                "{seg}/config/clades_genome.tsv",
        include =               "config/include.txt",
        meta =                  "data/metadata.tsv",
        meta_collab =           "data/meta_collab.tsv", ###TODO: Add an empty tsv file to this path or metadata for one of your sequences
        meta_manually_added = "data/meta_manually_added.tsv" # Published metadata


files = rules.files.input

##############################
# Download from NBCI Virus with ingest snakefile
###############################

rule fetch:  ####TODO: go to ingest readme and read through instructions and files needed!
    input:
        dir = "ingest"
    output:
        sequences="data/sequences.fasta",
        metadata=files.meta
    params:
        seq="ingest/data/sequences.fasta",
        meta="ingest/data/metadata.tsv"
    shell:
        """
        cd {input.dir} 
        snakemake --cores 9 all
        cd ../
        cp -u {params.seq} {output.sequences}
        cp -u {params.meta} {output.metadata}
        """

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
        metadata = files.meta,  # Path to input metadata file
        meta_collab = files.meta_collab,  # Data shared with us by collaborators
        meta_manually_added = files.meta_manually_added
    params:
        strain_id_field=config["id_field"],
        date_fields=config["curate"]["date_fields"],
        expected_date_formats=config["curate"]["expected_date_formats"],
    output:
        metadata = "data/curated/meta_public.tsv",  # Final output file for NCBI metadata
        meta_collab="data/curated/meta_collab.tsv",  # Curated collaborator metadata
        final_metadata="data/final_meta.tsv"  # Final merged output file
    shell:
        """
        # Normalize strings for publication metadata
        augur curate normalize-strings \
            --id-column {params.strain_id_field} \
            --metadata {input.metadata} \
        | augur curate format-dates \
            --date-fields {params.date_fields} \
            --no-mask-failure \
            --expected-date-formats {params.expected_date_formats} \
            --id-column {params.strain_id_field} \
            --output-metadata {output.metadata}
        
        # Normalize strings and format dates for collab metadata
        augur curate normalize-strings \
            --id-column {params.strain_id_field} \
            --metadata {input.meta_collab} \
        | augur curate format-dates \
            --date-fields {params.date_fields} \
            --no-mask-failure \
            --expected-date-formats {params.expected_date_formats} \
            --id-column {params.strain_id_field} \
            --output-metadata {output.meta_collab}
        
        # Merge curated metadata
        augur merge --metadata metadata={output.metadata} meta_collab={output.meta_collab} meta_manually_added={input.meta_manually_added} \
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata {output.final_metadata}
        """

##############################
# Update strain names
# If strain name == accession -> fetching real strain names from genbank
# Depending on how many sequences you have, it will run for a long time! >30min. Comment out to skip!
###############################

rule update_strain_names:
    message:
        """
        Updating strain information in metadata.
        """
    input:
        file_in =  rules.curate.output.final_metadata
    params:
        backup = "data/strain_names_previous_run.tsv" ####TODO: provide an empty file for first run
    output:
        file_out = "data/updated_strain_names.tsv"
    shell:
        """
        time bash scripts/update_strain.sh {input.file_in} {params.backup} {output.file_out}
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
        seqs_to_blast = rules.fetch.output.sequences
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
        input_seqs = rules.fetch.output.sequences
    output:
        sequences = "{seg}/results/sequences.fasta"
        
    params:
        range = "{seg}",  # Determines which protein (or whole genome) is processed
        # min_length = lambda wildcards: {"protein_xy": 2400, "protein_ab": 780, "protein_cd": 180, "protein_ef": 450, "genome": 6000}[wildcards.seg],  # Min length ## DONE replace whole_genome with genome, min length = 6000
        # max_length = lambda wildcards: {"protein_xy": 4000, "protein_ab": 1300, "protein_cd": 300, "protein_ef": 750, "genome": 8000}[wildcards.seg]  # Max length ## DONE replace whole_genome with genome, max length = 8000
        min_length = lambda wildcards: {"vp1": 700, "whole_genome": 6000}[wildcards.seg],  # Min length ## DONE replace whole_genome with genome, min length = 6000
        max_length = lambda wildcards: {"vp1": 900, "whole_genome": 8000}[wildcards.seg]    
    shell:
        """
        python scripts/blast_sort.py --blast {input.blast_result} \
            --seqs {input.input_seqs} \
            --out_seqs {output.sequences} \
            --range {params.range} \
            --min_length {params.min_length} \
            --max_length {params.max_length}

        # rm -r temp
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
        metadata =rules.curate.output.final_metadata,
        exclude = files.dropped_strains,
        include = files.include
    output:
        sequences = "{seg}/results/filtered.fasta"
    params:
        group_by = "country",
        sequences_per_group = 4000, # add a limit per group
        strain_id_field= config["id_field"],
        min_date = 1980,  # add a reasonable min date
        min_length = lambda wildcards: {"vp1": 700, "whole_genome": 6000}[wildcards.seg], 
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
            --output {output.sequences} \
            --min-length {params.min_length} 
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
        Aligning sequences to {input.reference} using Nextalign.
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        alignment = "{seg}/results/aligned.fasta"

    params:
            nuc_mismatch_all = 10,
            nuc_seed_length = 30
    shell:
        """
        nextclade run \
        {input.sequences}  \
        --input-ref {input.reference}\
        --allowed-mismatches {params.nuc_mismatch_all} \
        --min-length {params.nuc_seed_length} \
        --include-reference false \
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

    threads: 9
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
        metadata = rules.curate.output.final_metadata,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        tree = "{seg}/results/tree.nwk",
        node_data = "{seg}/results/branch_lengths.json"
    params:
        coalescent = "opt",
        rooting = lambda wildcards: {"vp1": "mid_point", "whole_genome": "mid_point"}[wildcards.seg], # "FR796483 GQ214176"
        date_inference = "marginal",
        clock_filter_iqd = 3, # set to 6 if you want more control over outliers
        strain_id_field = config["id_field"],
        clock_rate = 0.004, # remove for estimation by augur; check literature
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
        metadata =rules.curate.output.final_metadata
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
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
        """

#########################
#  EXPORT
#########################
rule export:
    message: "Creating auspice JSONs"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.curate.output.final_metadata,
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
        coloring_fields = "region country published_clade"
    output:
        auspice_json = "auspice/cva10_{seg}.json"
        
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} \
                {input.aa_muts} {input.clades} \
            --colors {input.colors} \
            --color-by-metadata {params.coloring_fields} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """
        
rule clean:
    message: "Removing directories: {params}"
    params:
        "*/results/*",
        "auspice/*",
        "temp/*"
    shell:
        "rm -rfv {params}"
