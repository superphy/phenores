
import os

#################################################################################
# FUNCTIONS                                                                     #
#################################################################################

def OPJ(*args):
    path = os.path.join(*args)
    return os.path.normpath(path)


#################################################################################
# GLOBALS                                                                       #
#################################################################################


#################################################################################
# RULES                                                                         #
#################################################################################


rule all:
    input:
        "data/interim/streptomycin_population_groups.csv"


rule str_files:
    input:
        "data/raw/Updated_GenotypicAMR_Master.xlsx"
    output:
        "data/interim/streptomycin_fasta_files.txt"
    params:
        fastadir="data/raw/genomes"
    run:
        import pandas as pd
        import numpy as np

        amrdf = pd.read_excel(input[0])

        amrdf = amrdf.replace(r'\s+', np.nan, regex=True)
        amrdf = amrdf.replace(r'-', np.nan, regex=True)

        # Drop Ecoli & friends
        amrdf = amrdf[amrdf['spcecies'] == 'enterica',:]
        sdf = amrdf[['run', 'phenotypic_streptomycin']].dropna()

        with open(output[0], 'w') as outfh:
            for index, row in sdf.iterrows():
                filepath = '{}/{}.fasta'.format(params.fastadir,row['run'])
                if not os.path.exists(filepath):
                    raise Exception('File does not exist: {}'.format(filepath))
                outfh.write(filepath+"\n")


rule str_mash:
    # Compute mash distance between all genomes
    input:
        "data/interim/streptomycin_fasta_files.txt"
    output:
        "data/interim/streptomycin_mash_distances.txt"
    params:
        sketchfile="data/interim/streptomycin_sketch",
        k=21, # K-mer length
        p=12, # Threads
        S=42, # Seed
        s=10000 # Sketch size
    shell:
        """
        mash sketch -k {params.k} -p {params.p} -s {params.s} -S {params.S} -l {input} -o {params.sketchfile} && \
        mash dist -p {params.p} {params.sketchfile}.msh {params.sketchfile}.msh > {output}
        """

rule str_dist:
    # Convert mash distance output into matrix
    input:
        "data/interim/streptomycin_mash_distances.txt"
    output:
        "data/interim/streptomycin_distance_matrix.csv"
    run:
        import pandas as pd
        import os

        ddf = pd.read_csv(input[0], sep='\t', header=None, names=['f1','f2','d','p','h'])

        ddf['s1'] = ddf.f1.apply(lambda x: os.path.splitext(os.path.os.path.basename(x))[0])
        ddf['s2'] = ddf.f2.apply(lambda x: os.path.splitext(os.path.os.path.basename(x))[0])
        
        dist = ddf.pivot(index='s1', columns='s2', values='d')
        dist.to_csv(output[0])


rule str_groups:
    # Assign genomes to clusters based on mash distances
    input:
        "data/interim/streptomycin_distance_matrix.csv",
        "data/raw/Updated_GenotypicAMR_Master.xlsx"
    output:
        "data/interim/streptomycin_population_groups.csv"
    params:
        k=10 # Number of clusters
    script: "src/data/make_mash_groups.R"






