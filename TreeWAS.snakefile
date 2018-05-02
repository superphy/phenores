
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
DATASETS='gene genome'.split()
#DATASETS='gene genome kmer dsm'.split()
METADATAFILE="data/interim/streptomycin_population_groups.csv"
METADATAFILE="data/interim/test_pop.csv"


#################################################################################
# RULES                                                                         #
#################################################################################

rule all:
    input:
        expand("data/interim/treewas/{fold}/{ds}/treewas_results.rdata", fold=FOLDS, ds=DATASETS)


rule binary:
    # Convert roary/piggy outputs to binary matrix
    input:
        METADATAFILE,
        "data/interim/roary/gene_presence_absence.csv",
        "data/interim/roary/IGR_presence_absence.csv",
    output:
        "data/interim/roary/gene_presence_absence_matrix.csv",
        "data/interim/roary/gene_and_igr_presence_absence_matrix.csv",
    run:
        import pandas as pd

        # Rtabs in Piggy are missing cluster names, so we will do it the hard way
        genedf = pd.read_csv(input[1], sep=',', header=0, index_col=0, na_values="", dtype=str)
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


rule filter_genes:
    # Filter original gene_presence_absence_matrix.csv to include only relevent training genomes
    input:
        METADATAFILE,
        "data/interim/roary/gene_presence_absence_matrix.csv"
    output:
        expand("data/interim/treewas/{fold}/gene/features.csv", fold=FOLDS)
    run:
        from Bio import SeqIO
        from contextlib import ExitStack
        import pandas as pd
        import numpy as np
        import os

        # Load groups
        df = pd.read_csv(input[0], sep=',', header=0, index_col=0)

        trainingset = {}
        for r in df.itertuples():
            g = r.sample
            ts = r.fold//2
            trainingset[g] = ts

        # Open files for writing
        with ExitStack() as stack:
            files = [stack.enter_context(open(fname, 'w')) for fname in output]

            with open(input[1], 'r') as infh:
                # Add headers
                header = infh.readline()
                for fh in files:
                    fh.write(header)

                # Add rest to correct file
                for line in infh:
                    sample = line.split(',', 1)[0]
                    if sample in trainingset:
                        ds = trainingset[sample]
                        for i in FOLDS:
                            if i != ds:
                                fh = files[i]
                                fh.write(line)


rule filter_genomes:
    # Filter original gene_presence_absence_matrix.csv to include only relevent training genomes
    input:
        METADATAFILE,
        "data/interim/roary/gene_and_igr_presence_absence_matrix.csv"
    output:
        expand("data/interim/treewas/{fold}/genome/features.csv", fold=FOLDS)
    run:
        from Bio import SeqIO
        from contextlib import ExitStack
        import pandas as pd
        import numpy as np
        import os

        # Load groups
        df = pd.read_csv(input[0], sep=',', header=0, index_col=0)

        trainingset = {}
        for r in df.itertuples():
            g = r.sample
            ts = r.fold//2
            trainingset[g] = ts

        # Open files for writing
        with ExitStack() as stack:
            files = [stack.enter_context(open(fname, 'w')) for fname in output]

            with open(input[1], 'r') as infh:
                # Add headers
                header = infh.readline()
                for fh in files:
                    fh.write(header)

                # Add rest to correct file
                for line in infh:
                    sample = line.split(',', 1)[0]
                    if sample in trainingset:
                        ds = trainingset[sample]
                        for i in FOLDS:
                            if i != ds:
                                fh = files[i]
                                fh.write(line)


rule filter_pheno:
    # Filter original resistance phenotype data to include only relevent training genomes
    input:
        METADATAFILE,
    output:
        expand("data/interim/treewas/{fold}/phenotypes.csv", fold=FOLDS)
    run:
        from Bio import SeqIO
        from contextlib import ExitStack
        import pandas as pd
        import numpy as np
        import os

        # Load groups
        df = pd.read_csv(input[0], sep=',', header=0, index_col=0)

        trainingset = {}
        pheno = {}
        for r in df.itertuples():
            g = r.sample
            ts = r.fold//2
            trainingset[g] = ts
            pheno[g] = r.resistant

        # Open files for writing
        with ExitStack() as stack:
            files = [stack.enter_context(open(fname, 'w')) for fname in output]

            for fh in files:
                fh.write("sample,resistant\n")

            for g in trainingset:
                ds = int(trainingset[g])
                for i in FOLDS:
                    if i != ds:
                        fh = files[i]
                        fh.write("{},{}\n".format(g, pheno[g]))        


rule filter_aln:
    # Filter original core_gene_alignment.aln to include only relevent training genomes
    input: 
        METADATAFILE,
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


rule tree:
    # Build phylogenetic tree 
    input:
       "data/interim/treewas/{fold}/train_core.aln"
    output:
        "data/interim/treewas/{fold}/tree.nwk"
    shell:
        "FastTree -gtr -nt -nosupport -fastest -noml {input} > {output}"


# rule tree:
#     # Build phylogenetic tree 
#     input:
#        "data/interim/treewas/{fold}/train_core.aln"
#     output:
#         "data/interim/treewas/{fold}/tree.nwk"
#     shell:
#         "clearcut -k --DNA --alignment --in={input} --out={output}"


rule treewas_genes:
    # Run treewas with gene presence / absence run
    input:
        "data/interim/treewas/{fold}/gene/features.csv",
        "data/interim/treewas/{fold}/phenotypes.csv",
        "data/interim/treewas/{fold}/tree.nwk"
    output:
        "data/interim/treewas/{fold}/gene/treewas_results.rdata",
        "data/interim/treewas/{fold}/gene/treewas_plots.pdf"
    script: 
        "src/gwas/run_treewas.R"


rule treewas_genomes:
    # Run treewas with gene presence / absence run
    input:
        "data/interim/treewas/{fold}/genome/features.csv",
        "data/interim/treewas/{fold}/phenotypes.csv",
        "data/interim/treewas/{fold}/tree.nwk"
    output:
        "data/interim/treewas/{fold}/genome/treewas_results.rdata",
        "data/interim/treewas/{fold}/genome/treewas_plots.pdf"
    script: 
        "src/gwas/run_treewas.R"


        
                





