#!/usr/bin/env python

"""make_kmer_groups.py

Using kmer distribution, compute cluster groups

"""

from kmer import KmerTable

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2018, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


def main():

    kmertable = KmerTable(snakemake.input['dir'], None)
    kmertable.sample(0.01)
    

if __name__ == '__main__':
    
    main()
