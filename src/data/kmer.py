#!/usr/bin/env python

"""kmer.py

Classes to convert dsk counts to sklearn objects

"""

import logging
import os
import numpy as np
import tables as pt
from subprocess import call
from scipy.sparse import csr_matrix

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2018, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


class COO(pt.IsDescription):
    row       = pt.UInt32Col()
    col       = pt.UInt32Col()

class Kmers(pt.IsDescription):
    kmer      = pt.StringCol(32)
    idx       = pt.UInt32Col()


class KmerTable:
    """KmerTable Class

    Pytable Kmer storage class

    """

    def __init__(self, hdf5_filepath, dsk_opts=None):

        self.file = hdf5_filepath
        dsk_defaults = {
            '-kmer-size': 31,
            '-nb-core': 1,
            '-out-tmp': '/tmp',
            '-out-dir': '/tmp',
            '-verbose': 1,
            '-abundance-min': 1
        }

        self.dsk_opts = dsk_opts or dsk_defaults
        self.dsk_cmd = 'dsk'
        self.dsk2ascii_cmd = 'dsk2ascii'
        self._tmp_kmer_counter = 0
        self._tmp_genome_counter = 0
       
    
    def populate(self, fasta_file_list):
        """Create pytable objects and fill with kmer data

        Overwrites any existing files

        Args:
            fasta_file_list(list): List of fasta file paths, one per line

        Returns:
            None

        """

        h5file = pt.open_file(self.file, mode='w')

        self.init(h5file)
        
        with open(fasta_file_list, 'r') as infh:
            for f in infh:
                f = f.strip()
                self._store(f, h5file)

        h5file.close()


    def init(self, h5file):
        """Create pytable objects

        Args:
            h5file(pytables.File)

        Returns:
            None

        """

        h5file.create_table(h5file.root, 'indexes', COO, "COO Indexes", expectedrows=1000000000)
        h5file.create_table(h5file.root, 'kmers', Kmers, "Kmer Lookup Table", expectedrows=10000000)
        genomeid_atom = pt.StringAtom(itemsize=11)
        h5file.create_earray(h5file.root, 'genomes', genomeid_atom, (0,), "Genome order", expectedrows=1500)
        
        h5file.flush()



    def add(self, fasta_file):
        pass


    def _store(self, fasta_file, h5file):
        """Run dsk and save to pytable file

        Args:
            fasta_file(str): Fasta file containing genome to add to table

        Returns:
            None

        """
        genomeid = os.path.splitext(os.path.basename(fasta_file))[0]
        dsk_file = os.path.join(
            self.dsk_opts['-out-tmp'],
            genomeid + '_dsk.h5'
        )
        kmer_file = os.path.join(
            self.dsk_opts['-out-tmp'],
            genomeid + '_dsk_kmers.txt'
        )

        # Run dsk
        self._rundsk(fasta_file, dsk_file, kmer_file)

        # Load into pytables
        self._loaddsk(kmer_file, h5file, genomeid)

        # Remove tmp files
        os.remove(dsk_file)
        #os.remove(kmer_file)


    def _rundsk(self, fasta_file, dsk_file, kmer_file):
        cmd = [self.dsk_cmd]
        for key,val in self.dsk_opts.items():
            cmd.extend([str(key), str(val)])
        cmd.extend(['-file', fasta_file, '-out', dsk_file])
        print(' '.join(cmd))
        call(cmd)


        cmd2 = [self.dsk2ascii_cmd]
        cores = str(self.dsk_opts['-nb-cores'])
        cmd2.extend(['-nb-cores', cores, '-file', dsk_file, '-out', kmer_file])
        print(' '.join(cmd2))
        call(cmd2)


    def _loaddsk(self, kmer_file, h5file, genomeid):

        h5file.root.genomes.append(np.array([genomeid], dtype='U11'))
        genome_idx = self._tmp_genome_counter
        self._tmp_genome_counter += 1
        

        kmertable = h5file.root.kmers
        indextable = h5file.root.indexes
        
        with open(kmer_file, 'r') as infh:
            for l in infh:
                (kmerstr, freq) = l.split()
                result = [row['idx'] for row in kmertable.where('kmer == kmerstr')]

                kmer_idx = self._tmp_kmer_counter
                if not result or len(result) == 0:
                    # New kmer
                    kmerrow = kmertable.row
                    kmerrow['kmer'] = kmerstr
                    kmerrow['idx'] = kmer_idx
                    kmerrow.append()
                    self._tmp_kmer_counter += 1

                elif len(result) > 1:
                    raise Exception('Multiple identical kmers')

                else:
                    # Existing kmer
                    kmer_idx = result[0]

                # Save new data point
                indexrow = indextable.row
                indexrow['row'] = genome_idx
                indexrow['col'] = kmer_idx
                indexrow.append()

        h5file.flush()

        
    def csr_matrix(self):
        """Return scipy csr_matrix representation of kmer table

        Args:
            None

        Returns tuple:
            [1]: csr_matrix
            [2]: numpy array with genome/row index
            [3]: numpy array with kmer/col index

        """

        h5file = pt.open_file(self.file, mode='r')

        n = h5file.root.indexes.nrows

        cm = csr_matrix(([1]*n, (h5file.root.indexes.cols.row, h5file.root.indexes.cols.col)))

        return (cm, h5file.root.genomes.read(), h5file.root.kmers.read(field='kmer'))
        


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.DEBUG, format=log_fmt)


    #main(snakemake.input[0])
