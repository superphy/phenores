
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
        "data/interim/streptomycin/kmers/kmer_matrix.txt"


rule streptokmers:
    input:
        "data/interim/streptomycin_fasta_files.txt"
    output:
        dir="data/interim/streptomycin/kmers/kmer_ohe_zarr"
    params:
        k=11,
        mins=1,
        v=0,
        cores=12,
        d="data/interim/dsk"
    script: 
        "src/data/make_kmer_table.py"


rule filter:
    input:
        "data/interim/streptomycin/kmers/kmer_ohe_zarr"
    output:
        "data/interim/streptomycin/kmers/kmer_matrix.txt"
    run:
        import zarr
        import pandas as pd
        import dask.array as da
        import dask.dataframe as dd
        import numpy as np

        store = zarr.DirectoryStore(input[0])
        group = zarr.hierarchy.group(store=store,overwrite=False)
        za = group['ohe']

        a = da.from_array(za, chunks=za.chunks)

        with dask.set_options(scheduler='threads'):

            # Filter by counts
            frequencies = a.sum(axis=0)
            mask = (frequencies > 14) & (frequencies < 1392)
            a = a[:,mask]

            df = dd.from_array(a.compute())
            df.to_csv()
            


        













