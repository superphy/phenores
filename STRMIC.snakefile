
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
        "data/interim/str_mic_class_dataframe.pkl", 
        "data/interim/str_mic_class_order_dict.pkl"

rule mics:
    input: 
        "data/raw/GenotypicAMR_Master.csv"
    output:
        "data/interim/str_mic_class_dataframe.pkl", "data/interim/str_mic_class_order_dict.pkl"
    script:
        "src/data/bin_str_mics.py"


