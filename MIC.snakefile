
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

rule all:
    input:
        "data/interim/mic_class_dataframe2.pkl", 
        "data/interim/mic_class_order_dict2.pkl"

rule mics:
    input: 
        "data/raw/Updated_GenotypicAMR_Master.xlsx"
    output:
        "data/interim/mic_class_dataframe2.pkl", 
        "data/interim/mic_class_order_dict2.pkl"
    script:
        "src/data/bin_mics.py"


