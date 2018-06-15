import dotenv
import os
import glob

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
PROJECT_DIR = './'

# dotenv project variables
dotenv_path = OPJ(PROJECT_DIR, ".env")
if not os.path.exists(dotenv_path):
    raise Exception("Missing .env file: {}".format(dotenv_path))
dotenv.load_dotenv(dotenv_path)

FOLDS=range(5)
#DATASETS='strings kmers'.split()
DATASETS='strings'.split()
METADATAFILE="data/interim/streptomycin_population_groups.csv"
#METADATAFILE="data/interim/test_pop.csv"

##################################################################################
## TARGETS                                                                       #
##################################################################################
# import pandas as pd
# GENOMES=[]
# df = pd.read_csv("data/interim/streptomycin_population_groups.csv", sep=',', header=0, index_col=0)
# for r in df.itertuples():
#     GENOMES.append(r.sample)

#GENOMES=GENOMES[0:9]

METADATAFILE="data/interim/streptomycin_population_groups.csv"

#################################################################################
# RULES                                                                         #
#################################################################################

#localrules: all


rule all:
    input:
        expand("data/interim/seer/{fold}/{ds}/seer_results.txt", fold=FOLDS, ds=DATASETS)


rule fsmlite:
	# Takes a long time
	input:
		"data/interim/streptomycin_fasta_files.txt"
	params:
		tmp="data/interim/fsmlite/tmp/fsml_index",
		smin=14,
		smaj=1392
	output:
		"data/interim/fsmlite/streptomycin_kmers.txt"
	shell:
		"../fsm-lite/fsm-lite -v -l {input[0]} -t {params.tmp} -s {params.smin} -S {params.smaj} > {output[0]}"


rule gzip:
	input:
		"data/interim/fsmlite/streptomycin_kmers.txt"
	output:
		"data/interim/fsmlite/streptomycin_kmers.txt.gzip"
	shell:
		"""
		gzip {input[0]}
		"""





