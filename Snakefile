
import os

configfile: "config.yaml"

#################################################################################
# FUNCTIONS                                                                     #
#################################################################################

def OPJ(*args):
    path = os.path.join(*args)
    return os.path.normpath(path)


#################################################################################
# GLOBALS                                                                       #
################################################################################# PROJECT_NAME = 'phenores' PROJECT_DIR = OPJ(os.path.dirname(__file__), os.pardir)

#################################################################################
# RULES                                                                         #
#################################################################################


# rule roary:
# 	input:
# 		"data/external/AAC.02140-16_zac003175944sd1.csv"
# 	output:
# 		"data/interim/accession_map.txt"
# 	shell:
# 		"""
# 		cut -f2 {input} | sort > data/interim/PMC5328538_original_genome_accessions.sort
# 		grep -v -P "^SRR"  data/interim/PMC5328538_original_genome_accessions.sort > data/interim/PMC5328538_assembly_ids.txt
# 		cut -f3,5 data/external/PRJNA242614_AssemblyDetails.txt | sort | join - data/interim/PMC5328538_assembly_ids.txt > data/interim/PMC5328538_assembly_biosample_ids.txt
# 		(
# 		cd data/interim
# 		./get_sra_accession.sh
# 		)
# 		(
# 		cd data/interim
# 		./merge.sh
# 		)
# 		sort -t' ' -k3 data/interim/PMC5328538_assembly_biosample_sra_ids.txt | join -t' ' -11 -23 -a1 -o1.1,2.1 data/interim/PMC5328538_sra_ids.txt - > {output}
# 		


rule mics:
	input:
		"data/raw/GenotypicAMR.xlsx"
	script:
		"src/data/bin_mics.py"

rule countkmers:
    input:
        "data/raw/genomes/"
    output:
        expand("data/interim/db.k{k}.l{l}", k=config["k"], l=config["l"])
    run:
        from kmerprediction import kmer_counter
        all_files = []
        for f in input:
            files = [f + x for x in os.listdir(f)]
            all_files.extend(files)
        kmer_counter.count_kmers(config["k"], config["l"], files, output, True)

rule gatherdata:
    input:
        "data/raw/genomes/",
        "data/raw/GenotypicAMR.xlsx",
        expand("data/interim/db.k{k}.l{l}", k=config["k"], l=config["l"]),
    output:
        expand("data/interim/train_test_split.k{k}.l{l}/train_data.pkl", k=config["k"], l=config["l"]),
        expand("data/interim/train_test_split.k{k}.l{l}/train_labels.pkl", k=config["k"], l = config["l"]),
        expand("data/interim/train_test_split.k{k}.l{l}/test_data.pkl", k=config["k"], l=config["l"]),
        expand("data/interim/train_test_split.k{k}.l{l}/test_labels.pkl", k=config["k"], l=config["l"]),
        expand("data/interim/train_test_split.k{k}.l{l}/label_encoder.pkl", k=config["k"], l=config["l"]),
    script:
        "src/data/gatherdata.py"

rule preprocessdata:
    input:
        expand("data/interim/train_test_split.k{k}.l{l}/train_data.pkl", k=config["k"], l=config["l"]),
        expand("data/interim/train_test_split.k{k}.l{l}/train_labels.pkl", k=config["k"], l = config["l"]),
        expand("data/interim/train_test_split.k{k}.l{l}/test_data.pkl", k=config["k"], l=config["l"]),
        expand("data/interim/train_test_split.k{k}.l{l}/test_labels.pkl", k=config["k"], l=config["l"])
    output:
        "data/processed/train_data.pkl",
        "data/processed/train_labels.pkl",
        "data/processed/test_data.pkl",
        "data/processed/test_labels.pkl",
    script:
        "src/data/preprocessdata.py"

rule trainneuralnet:
    input:
        "data/processed/train_data.pkl",
        "data/processed/train_labels.pkl",
        "data/processed/test_labels.pkl"
    output:
        "models/neural_net.h5",
    script:
        "src/models/train_neural_net.py"

rule testneuralnet:
    input:
        "models/neural_net.h5",
        "data/processed/test_data.pkl",
        "data/processed/test_labels.pkl"
    script:
        "src/models/test_neural_net.py"


