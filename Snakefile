
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
# 		"""


rule mics:
	input: 
		"data/raw/GenotypicAMR.csv"
	output:
		"data/interim/mic_class_dataframe.pkl", "data/interim/mic_class_order_dict.pkl"
	script:
		"src/data/bin_mics.py"


rule streptofiles:
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
		"data/interim/streptomycin_kmers.h5"
	params:
		k=31,
		mins=1,
		d="data/interim/dsk",
		cores=8,
		v=1
	script: "src/data/make_kmer_table.py"
	

rule test_streptokmers:
