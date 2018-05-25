#!/usr/bin/env python

"""bin_str_mics.py

Convert Streptomycin MIC values into distinct classes. Save processed
classes as pandas dataframes.

"""

import os
import logging
import re
import numpy as np
import pandas as pd

from dotenv import find_dotenv, load_dotenv
from sklearn.externals import joblib

from mic import MICPanel, MGPL

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2018, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"

# Manually define MIC ranges due to mixing of different systems
mic_ranges = {
    'phenotypic_streptomycin.1': {
        'top': '>=64.0000',
        'bottom': '<=2.0000',
    },
}


def main(excel_filepath):
    """ Runs data processing scripts to turn MIC text values from Excel input data
        into class categories

    Args:
        excel_filepath: Metadata Excel file. MIC columns will have prefix 'MIC_'

    """
    logger = logging.getLogger(__name__)
    logger.info('Streptomycin MIC binning')

    drugs = ["phenotypic_streptomycin.1"]
    convert = { key: lambda x: str(x) for key in drugs }
   
    micsdf = pd.read_csv(excel_filepath, sep='\t', usecols=['run'] + drugs, skipfooter=1, skip_blank_lines=False, converters=convert)
    micsdf = micsdf.set_index('run')
    
    classes = {}
    class_orders = {}
    for col in micsdf:
        logger.debug('Creating MIC panel for {}'.format(col))
        class_labels, class_order = bin(micsdf[col], col)
        match = re.search(r'^phenotypic_(\w+)(?:\.1)?$', col)
        drug = match.group(1)
        classes[drug] = pd.Series(class_labels, index=micsdf.index)
        class_orders[drug] = class_order

        logger.debug("Final MIC distribution:\n{}".format(classes[drug].value_counts()))

    c = pd.DataFrame(classes)

    data_dir = os.environ.get('PRDATA')
    cfile = os.path.join(data_dir, 'interim', 'str_mic_class_dataframe.pkl')
    cofile = os.path.join(data_dir, 'interim', 'str_mic_class_order_dict.pkl')
    joblib.dump(c, cfile)
    joblib.dump(class_orders, cofile)



def bin(mic_series, drug):

    logger = logging.getLogger(__name__)

    panel = MICPanel()
    values, counts = panel.build(mic_series)

    logger.debug('Panel value frequency:')
    for m,c in zip(values, counts):
        logger.debug("{}, {}".format(m, c))


    # Initialize MIC bins
    panel.set_range(mic_ranges[drug]['top'], mic_ranges[drug]['bottom'])

    # Iterate through MIC values and assign class labels
    logger.debug('MIC values will be mapped to: {}'.format(panel.class_labels))
    classes = []
    for m in mic_series:
        mgpl = MGPL(m)
        mlabel = str(mgpl)

        if mgpl.isna:
            classes.append(np.nan)
        else:
            if not mlabel in panel.class_mapping:
                raise Exception('Mapping error')
            else:
                classes.append(panel.class_mapping[mlabel])


    return(classes, panel.class_labels)




if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.DEBUG, format=log_fmt)

    # Load environment
    project_dir = os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)
    load_dotenv(find_dotenv())

    main(snakemake.input[0])
