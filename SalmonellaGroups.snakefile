
import os

#################################################################################
# FUNCTIONS                                                                     #
#################################################################################


#################################################################################
# GLOBALS                                                                       #
#################################################################################


#################################################################################
# RULES                                                                         #
#################################################################################


rule all:
    input:
        "data/interim/mash_population_groups.csv"


rule mash:
    # Compute mash distance between all genomes
    input:
        "data/interim/fasta_files_list_no_ecoli.txt"
    output:
        "data/interim/mash_distances.txt"
    params:
        sketchfile="data/interim/mash_sketch",
        k=21, # K-mer length
        p=12, # Threads
        S=42, # Seed
        s=10000 # Sketch size
    shell:
        """
        mash sketch -k {params.k} -p {params.p} -s {params.s} -S {params.S} -l {input} -o {params.sketchfile} && \
        mash dist -p {params.p} {params.sketchfile}.msh {params.sketchfile}.msh > {output}
        """

rule dist:
    # Convert mash distance output into matrix
    input:
        "data/interim/mash_distances.txt"
    output:
        "data/interim/mash_distance_matrix.csv"
    run:
        import pandas as pd
        import os

        ddf = pd.read_csv(input[0], sep='\t', header=None, names=['f1','f2','d','p','h'])

        ddf['s1'] = ddf.f1.apply(lambda x: os.path.splitext(os.path.os.path.basename(x))[0])
        ddf['s2'] = ddf.f2.apply(lambda x: os.path.splitext(os.path.os.path.basename(x))[0])
        
        dist = ddf.pivot(index='s1', columns='s2', values='d')
        dist.to_csv(output[0])


rule groups:
    # Assign genomes to clusters based on mash distances
    input:
        "data/interim/mash_distance_matrix.csv",
    output:
        "data/interim/mash_population_groups.csv"
    params:
        k=9 # Number of clusters
    script: "src/data/make_mash_salmonella_groups.R"






