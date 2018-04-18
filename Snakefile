
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

PROJECT_NAME = 'phenores'
PROJECT_DIR = OPJ(os.path.dirname(__file__), os.pardir)
FOLDS=range(5)

#################################################################################
# RULES                                                                         #
#################################################################################


# rule roary:
#   input:
#       "data/external/AAC.02140-16_zac003175944sd1.csv"
#   output:
#       "data/interim/accession_map.txt"
#   shell:
#       """
#       cut -f2 {input} | sort > data/interim/PMC5328538_original_genome_accessions.sort
#       grep -v -P "^SRR"  data/interim/PMC5328538_original_genome_accessions.sort > data/interim/PMC5328538_assembly_ids.txt
#       cut -f3,5 data/external/PRJNA242614_AssemblyDetails.txt | sort | join - data/interim/PMC5328538_assembly_ids.txt > data/interim/PMC5328538_assembly_biosample_ids.txt
#       (
#       cd data/interim
#       ./get_sra_accession.sh
#       )
#       (
#       cd data/interim
#       ./merge.sh
#       )
#       sort -t' ' -k3 data/interim/PMC5328538_assembly_biosample_sra_ids.txt | join -t' ' -11 -23 -a1 -o1.1,2.1 data/interim/PMC5328538_sra_ids.txt - > {output}
#       """

rule all:
    input:
        expand("data/interim/treewas/{fold}/train_core.nwk", fold=FOLDS)


rule mics:
    input: 
        "data/raw/GenotypicAMR.csv"
    output:
        "data/interim/mic_class_dataframe.pkl", "data/interim/mic_class_order_dict.pkl"
    script:
        "src/data/bin_mics.py"


rule str_files:
    input:
        "data/raw/GenotypicAMR_Master.xlsx"
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

        sdf = amrdf[['run', 'phenotypic_streptomycin']].dropna()

        with open(output[0], 'w') as outfh:
            for index, row in sdf.iterrows():
                filepath = '{}/{}.fasta'.format(params.fastadir,row['run'])
                if not os.path.exists(filepath):
                    raise Exception('File does not exist: {}'.format(filepath))
                outfh.write(filepath+"\n")


rule streptokmers:
    input:
        "data/interim/streptomycin_fasta_files.txt"
    output:
        dir="data/interim/streptomycin/kmers/"
    params:
        k=31,
        mins=1,
        v=0,
        cores=12,
        d="data/interim/dsk"
    script: "src/data/make_kmer_table.py"


# rule streptoclusters:
#   input:
#       dir="data/interim/streptomycin/kmers/"
#   output:
#       "data/interim/streptomycin_clusters.pkl"
#   script: "src/data/make_kmer_groups.py"


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
    # Assign strings to clusters based on mash distances
    input:
        "data/interim/streptomycin_distance_matrix.csv",
        "data/raw/GenotypicAMR_Master.xlsx"
    output:
        "data/interim/streptomycin_population_groups.csv"
    params:
        k=10 # Number of clusters
    script: "src/data/make_mash_groups.R"


rule str_neptune_files:
    # Prepare neptune inputs
    # Neptune and snakemake are incompatible, need to run separately
    input:
        "data/interim/streptomycin_population_groups.csv"
    output:
        expand('data/interim/neptune/{fold}/{folder}/', fold=FOLDS, folder='inclusive exclusive validation'.split())
    run:
        import pandas as pd
        import os

        df = pd.read_csv(input[0], sep=',', header=0, index_col=0)
        
        # Make directory
        for d in output:
            if not os.path.exists(d):
                os.makedirs(d)

        # Symlink fasta files to each directory
        for f in FOLDS:
            inclusive_path=OPJ("data","interim","neptune",str(f),"inclusive")
            exclusive_path=OPJ("data","interim","neptune",str(f),"exclusive")
            validation_path=OPJ("data","interim","neptune",str(f),"validation")
            for r in df.itertuples():
                file = os.path.basename(r.fasta)
                srcfile = os.path.abspath(r.fasta)
                if r.fold >= f*2 and r.fold <= f*2+1:
                    # Validation
                    os.symlink(srcfile, OPJ(validation_path,file))
                else:
                    if r.resistant:
                        os.symlink(srcfile, OPJ(inclusive_path,file))
                    else:
                        os.symlink(srcfile, OPJ(exclusive_path,file))


rule str_filter_aln:
    # Filter original core_gene_alignment.aln to include only relevent training genomes
    input: 
        "data/interim/streptomycin_population_groups.csv",
        "data/interim/roary/core_gene_alignment.aln"
    output: 
        expand("data/interim/treewas/{fold}/train_core.aln", fold=FOLDS)
    run:
        from Bio import SeqIO
        import pandas as pd
        import numpy as np
        import os

        # Load groups
        df = pd.read_csv(input[0], sep=',', header=0, index_col=0)

        # Load sequences
        sequences = {}
        for seq_record in SeqIO.parse(input[1], "fasta"):
            sequences[seq_record.id] = seq_record.seq
        
        # Output to relevent folders
        for fold in FOLDS:

            # Make directory
            d = OPJ("data","interim","treewas",str(fold))
            if not os.path.exists(d):
                os.makedirs(d)

            # Output training fasta file
            trainfile=OPJ(d, 'train_core.aln')
            with open(trainfile, 'w') as trainfh:
                
                for r in df.itertuples():
                    g = r.sample
                    if r.fold < fold*2 or r.fold > fold*2+1:
                        # Training set
                        seq = sequences[g]
                        if not seq:
                            raise Exception("Missing sequence for genome {}".format(g))
                        trainfh.write(">{}\n{}\n".format(g, seq))

# TREEWAS



rule str_core_tree:
    # Build phylogenetic tree 
    input:
       "data/interim/treewas/{fold}/train_core.aln"
    input:
        "data/interim/treewas/{fold}/train_core.nwk"
    shell:
        "FastTree -gtr -nt -nosupport -fastest -noml {input} > {output}"


rule str_gene_binary_matrix:
    # Convert roary/piggy outputs to binary matrix
    input:
        "data/interim/streptomycin_population_groups.csv",
        "data/interim/roary/gene_presence_absence.csv",
        "data/interim/roary/IGR_presence_absence.csv",
    output:
        "data/interim/roary/gene_presence_absence_matrix.csv",
        "data/interim/roary/gene_and_igr_presence_absence_matrix.csv",
    run:
        import pandas as pd

        # Rtabs in Piggy are missing cluster names, so we will do it the hard way
        genedf = pd.read_csv(input[1], sep=',', header=0, index_col=0, na_values="", dtype=str)
        inputs = "../data/interim/roary/IGR_presence_absence.csv"
        igrdf = pd.read_csv(input[2], sep=',', header=0, index_col=0, na_values="", dtype=str)
        sampledf = pd.read_csv(input[0], sep=',', header=0, index_col=0)

        genemat = genedf[sampledf["sample"]]
        genemat = genemat.applymap(lambda x: 0 if pd.isna(x) else 1)
        genemat = genemat.T

        igrmat = igrdf[sampledf["sample"]]
        igrmat = igrmat.applymap(lambda x: 0 if pd.isna(x) else 1)
        igrmat = igrmat.T

        genomemat = pd.concat([genemat, igrmat], axis=1)

        genemat.to_csv(output[0])
        genomemat.to_csv(output[1])


# rule str_treewas_genes:
#     # Run treewas with gene presence / absence run


        
                





