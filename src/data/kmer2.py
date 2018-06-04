#!/usr/bin/env python

"""kmer.py

Classes to convert dsk counts to sklearn objects

"""

import contextlib
import logging
import os
import plyvel
import shutil
import time
import zarr
from subprocess import call

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2018, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


class DskRunner:
    """Runs Dsk on a single fasta file

    """

    def __init__(self, counter_opts, fasta_file):

        dsk_defaults = {
            '-kmer-size': 31,
            '-nb-core': 1,
            '-out-tmp': '/tmp',
            '-out-dir': '/tmp',
            '-verbose': 1,
            '-abundance-min': 1
        }

        self.dsk_cmd = 'dsk'
        self.dsk2ascii_cmd = 'dsk2ascii'
        self.dsk_opts = counter_opts or dsk_defaults

        genomeid = os.path.splitext(os.path.basename(fasta_file))[0]
        self.dsk_file = os.path.join(
            self.dsk_opts['-out-tmp'],
            genomeid + '_dsk.h5'
        )
        self.kmer_file = os.path.join(
            self.dsk_opts['-out-tmp'],
            genomeid + '_dsk_kmers.txt'
        )
        self.fasta_file = fasta_file


    def run(self):

        fasta_file = self.fasta_file

        cmd = [self.dsk_cmd]
        for key,val in self.dsk_opts.items():
            cmd.extend([str(key), str(val)])
        cmd.extend(['-file', fasta_file, '-out', self.dsk_file])
        print(' '.join(cmd))
        call(cmd)

        cmd2 = [self.dsk2ascii_cmd]
        cores = str(self.dsk_opts['-nb-cores'])
        cmd2.extend(['-nb-cores', cores, '-file', self.dsk_file, '-out', self.kmer_file])
        print(' '.join(cmd2))
        call(cmd2)

        return self.iterkmers()


    def iterkmers(self):

        with open(self.kmer_file, encoding='utf-8') as file:
            for line in file:
                (kmerstr, freq) = line.split()
                yield kmerstr.encode()


    def close(self):

        os.remove(self.dsk_file)
        os.remove(self.kmer_file)


class DskCounter:
    """Wrapper around external kmer counting program Dsk

    Spawns DskRunner objects that can be run within a contextlib.closing
        block and therefore have any output files cleaned properly

    """

    def __init__(self, counter_opts):

        dsk_defaults = {
            '-kmer-size': 31,
            '-nb-core': 1,
            '-out-tmp': '/tmp',
            '-out-dir': '/tmp',
            '-verbose': 1,
            '-abundance-min': 1
        }

        self.dsk_opts = counter_opts or dsk_defaults


    def runner(self, fasta_file):

        return DskRunner(self.dsk_opts, fasta_file)



class KmerTable:
    """Kmer storage class

    """

    def __init__(self, zarr_storage_dir, leveldb_storage_dir, counter_object):

        self.counter = counter_object

        self.leveldb_filepath = leveldb_storage_dir
        self.zarr_filepath = zarr_storage_dir

        self._tmp_kmer_counter = 0
        self._tmp_genome_counter = 0
        self.buffer = 1000000
        self.chunk = 1000
       
    
    def enumerate(self, fasta_file_list):
        """Create pytable objects and fill with kmer data

        Overwrites any existing files

        Args:
            fasta_file_list(list): List of fasta file paths, one per line

        Returns:
            None

        """

        with open(fasta_file_list, 'r') as infh:
            fasta_files = infh.read().splitlines()
    
        self.n = len(fasta_files)

        # Kmer index lookup
        if os.path.exists(self.leveldb_filepath):
            shutil.rmtree(self.leveldb_filepath)
        db = plyvel.DB(self.leveldb_filepath, create_if_missing=True)

        # Zarr chunked dataframe
        store = zarr.DirectoryStore(self.zarr_filepath)
        group=zarr.hierarchy.group(store=store,overwrite=True)
        za=group.zeros('ohe', shape=(self.n, self.buffer), dtype='u1', chunks=(self.n, self.chunk))
        
        for f in fasta_files:
            self.count(f, db, za)

        db.close()


    def count(self, fasta_file, db, za):
        """Run dsk and save to pytable file

        Args:
            fasta_file(str): Fasta file containing genome to add to table

        Returns:
            None

        """
        
        with contextlib.closing(self.counter.runner(fasta_file)) as runner:
    
            genomeid = os.path.splitext(os.path.basename(fasta_file))[0]

            # Run dsk
            iterkmers = runner.run()

            # Load into pytables
            self.load(iterkmers, db, za, genomeid)


    def load(self, iterkmers, db, za, genomeid):

        wb = db.write_batch()
        tmp_n = 1000000

        start_time = time.time()

        genome_idx = self._tmp_genome_counter
        self._tmp_genome_counter += 1
        
        i = 0
        colbuffer = []
        for kmerbstr in iterkmers:

            t0 = time.time()

            # Logistic note: Each kmer is unique in a given genome file
            result = db.get(kmerbstr)  
               
            t1 = time.time()

            if not result:
                # New kmer
                kmer_idx = self._tmp_kmer_counter
                self._tmp_kmer_counter += 1

                wb.put(kmerbstr, str(kmer_idx).encode())

                if self._tmp_kmer_counter >= za.shape[1]:
                    ncols = za.shape[1]+self.buffer
                    za.resize(self.n, ncols)

            else:
                kmer_idx = int(result)

            t2 = time.time()

            # Save new data point
            colbuffer.append(kmer_idx)

            if len(colbuffer) >= tmp_n:
                s = len(colbuffer)
                try:
                    za.set_coordinate_selection(([genome_idx]*s, colbuffer), [1]*s)
                except IndexError:
                    mx = za.shape[1]
                    print("shape: {},{}".format(za.shape[0],za.shape[1]))
                    print("genome: {}".format(genome_idx))
                    print("overages colbuffer:\n{}\n".format(', '.join([str(c) for c in colbuffer if c >= mx ])))
                    raise
                colbuffer = []

            t3 = time.time()

            i += 1
            if i % 100000 == 0:
                print("{} completed. Timing: query {}, insert kmer row {}, insert index row {}".format(i,
                    t1-t0, t2-t1, t3-t2))
                diff = int(t3 - start_time)
                minutes, seconds = diff // 60, diff % 60
                print('\tElapsed type ' + str(minutes) + ':' + str(seconds).zfill(2))

        wb.write()

        

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.DEBUG, format=log_fmt)


    #main(snakemake.input[0])
