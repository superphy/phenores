#!/usr/bin/env python

"""make_kmer_table.py

Run dsk and save in pytable

"""

from kmer import KmerTable

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2018, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


def main():

    dsk_opts = {
        '-kmer-size': snakemake.params.k,
        '-nb-cores': snakemake.params.cores,
        '-abundance-min': snakemake.params.mins,
        '-verbose': snakemake.params.v,
        '-out-dir': snakemake.params.d,
        '-out-tmp': snakemake.params.d
    }

    kmer_table = KmerTable(snakemake.output[0], dsk_opts)

    kmer_table.populate(snakemake.input[0])
    


if __name__ == '__main__':
    
    main()
