#!/usr/bin/env python

"""Data Loading / Parsing Function Library


"""

import numpy as np
import os
import pandas as pd
import re
from sklearn.externals import joblib

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@canada.ca"



def process_metadata(input, sra_accession_map_file):
	"""Load raw supplementary file CSV data into pandas dataframes

	Divides meta data into logical datatypes and saves each datatype
	as separate pandas dataframes using joblib.dump

	Args:
		input(str): input filepath
		sra_accession_map_file(str): filepath to file that maps Assembly IDs to SRA accession IDs

	Returns:
		Saves:
			mics: MIC numeric values
			micssigns: Sign on MIC (e.g. >=)
			meta: contains columns serotype, year, source,
			marker: presence of AMR resistance gene markers


	"""

	df = pd.read_csv(input, delimiter='\t')

	# Map all assembly IDs to SRA accessinos
	accmap = pd.read_csv(sra_accession_map_file, delimiter=' ')
	df["Nucleotide accession"]
	df = df.set_index('Nucleotide accession')

	# MICs
	tmpdf = df[['AMP','AMC', 'FOX', 'CRO', 'TIO', 'GEN', 'FIS', 'SXT', 'AZM', 'CHL', 'CIP', 'NAL', 'TET']]
	micssign = tmpdf.applymap(lambda x: x.split()[0] if isinstance(x, str) else '=' )
	mic = tmpdf.applymap(lambda x: float(x.split()[1]) if isinstance(x, str) else float(x) )

	joblib.dump(mic, os.path.join(os.environ['PRDATA'], 'processed', 'mic.pkl'))
	joblib.dump(micssign, os.path.join(os.environ['PRDATA'], 'processed','micssign.pkl'))

	# Meta
	meta = df[['YEAR', 'SOURCE', 'SEROTYPE']]
	meta = meta.applymap(lambda x: x.strip().lower() if isinstance(x, str) else x )
	meta['SEROTYPE'] = meta['SEROTYPE'].apply(lambda x: re.sub(r' var\s', ' var. ', x))
	meta['SEROTYPE'] = meta['SEROTYPE'].apply(lambda x: re.sub(r'\s*:\s*', ':', x))
	joblib.dump(meta, os.path.join(os.environ['PRDATA'], 'processed', 'meta.pkl'))

	# Markers
	marker = df[['aminoglycoside', 'beta-lactam', 'macrolide', 'phenicol', 'quinolone', 'sulphonamide', 
		'tetracycline', 'trimethoprim',]]
	joblib.dump(marker, os.path.join(os.environ['PRDATA'], 'processed', 'marker.pkl'))

	return(df)


def load_metadata():
	"""Load pickle files into dataframes

	Returns:
		list: 
			mics: MIC numeric values
			micssigns: Sign on MIC (e.g. >=)
			meta: contains columns genus, serotype, year, source,
			marker: presence of AMR resistance gene markers

	"""

	mic = joblib.load(os.path.join(os.environ['PRDATA'], 'processed', 'mic.pkl'))
	micsign = joblib.load(os.path.join(os.environ['PRDATA'], 'processed','micssign.pkl'))
	meta = joblib.load(os.path.join(os.environ['PRDATA'], 'processed', 'meta.pkl'))
	marker = joblib.load(os.path.join(os.environ['PRDATA'], 'processed', 'marker.pkl'))

	return(mic, micsign, meta, marker)





