
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
        "data/interim/mic_class_dataframe2.pkl", 
        "data/interim/mic_class_order_dict2.pkl"

rule mics:
    input: 
        "data/raw/GenotypicAMR_Master.xlsx"
    output:
        "data/interim/mic_class_dataframe2.pkl", 
        "data/interim/mic_class_order_dict2.pkl"
    script:
        "src/data/bin_mics.py"

rule test:
    input:
        "data/raw/ResolvedCIPARS_SRL (1).xlsx"
    output:
        "data/interim/mic_class_dataframe_test.pkl",
        "data/interim/mic_class_order_dict_test.pkl"
    params:
        class_labels="data/interim/test_amp_classes.yaml"
    script:
        "src/data/bin_mics.py"


