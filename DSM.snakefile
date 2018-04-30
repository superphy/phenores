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

##################################################################################
## TARGETS                                                                       #
##################################################################################
import pandas as pd
GENOMES=[]
df = pd.read_csv("data/interim/streptomycin_population_groups.csv", sep=',', header=0, index_col=0)
for r in df.itertuples():
    GENOMES.append(r.sample)

#GENOMES=GENOMES[0:9]

#################################################################################
# RULES                                                                         #
#################################################################################

localrules: all, dsm_symlink, dsm_sample_file, dsm_run 

#rule all:
#    input:
#        expand("server-output.{H1}{H2}.txt.gz", H1='A C G T'.split(), H2='A C G T'.split())


rule all:
    input:
        "data/interim/dsm/streptomycin/samples.txt", 
        expand("data/interim/dsm/streptomycin/{genome}.fasta.fmi", genome=GENOMES)


rule dsm_symlink:
    output:
        "data/interim/dsm/streptomycin/{genome}.fasta"
    run:
        genomefile = os.path.abspath("data/raw/genomes/"+wildcards.genome+".fasta")
        dest = "data/interim/dsm/streptomycin/"+wildcards.genome+".fasta"
        os.symlink(genomefile, dest)

 
rule dsm_sample_file:
    output:
       "data/interim/dsm/streptomycin/samples.txt"
    params:
       genomestr="\n".join(GENOMES)
    shell:
       "echo -e '{params.genomestr}' > {output}"


rule dsm_build:
    input:
        "data/interim/dsm/streptomycin/{genome}.fasta"
    output:
        "data/interim/dsm/streptomycin/{genome}.fasta.fmi"
    params:
        cmd=os.environ.get("DSM_FRAMEWORK_PATH")+"builder"
    shell:
        "{params.cmd} -v {input}"

# NEED TO RUN THIS OUTSIDE OF A SNAKEMAKE FILE
# IT is already a dispatching script.
#rule dsm_run:
#    input:
#        "data/interim/dsm/streptomycin/samples.txt",
#        expand("data/interim/dsm/streptomycin/{genome}.fasta.fmi", genome=GENOMES),
#    output:
#        expand("server-output.{H1}{H2}.txt.gz", H1='A C G T'.split(), H2='A C G T'.split())
#    run:
#        shell("src/gwas/dsm-server.sh {input[0]} tmp_dsm_config"),
#        shell("sleep 30"),
#        shell("src/gwas/dsm-client.sh {input[0]} tmp_dsm_config")
        
