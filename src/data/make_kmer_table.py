#!/usr/bin/env python

"""make_kmer_table.py

Run dsk and save in pytable

"""

from kmer2 import KmerTable, DskCounter

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
        '-out-dir': snakemake.params.d,
        '-out-tmp': snakemake.params.d,
        '-verbose': snakemake.params.v,
    }

    dc = DskCounter(dsk_opts)
    kmertable = KmerTable(snakemake.output[0], snakemake.output[1], dc)
    kmertable.enumerate(snakemake.input[0]) 


if __name__ == '__main__':
    
    main()
