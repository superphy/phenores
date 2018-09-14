

#################################################################################
# FUNCTIONS                                                                     #
#################################################################################


#################################################################################
# GLOBALS                                                                       #
#################################################################################

def get_fasta(fl):

fasta_files = expand(shell("grep -v '^#' data/interim/fasta_files_list_no_ecoli.txt"))


#################################################################################
# RULES                                                                         #
#################################################################################


rule all:
    #input: "data/interim/gene_kmers/cmy2_kmer_ftests.tsv"
    input: "data/gene_kmers/nongrdi_blastdb.nin"


rule db:
    input:
        "data/interim/fasta_files_list_no_ecoli.txt"
    output:
        "data/gene_kmers/nongrdi_blastdb.nin"
    params:
        tmp="tmp.fasta",
        blastdb="data/gene_kmers/nongrdi_blastdb"
    shell:
        "cat $(grep -v '^#' {input[0]}) > {params.tmp}"
        "makeblastdb -in {params.tmp} -out {params.blastdb} -dbtype 'nucl' -max_file_sz '4GB'"

rule blast:
    input:
        "data/interim/fasta_files_list_no_ecoli.txt"
    output:
        "data/gene_kmers/cmy2_blast_results/log.txt"
    params:

    script:
        "src/data/make_kmer_table.py"