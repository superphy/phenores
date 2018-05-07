#!/usr/bin/env python

"""kmer.py

Classes to convert dsk counts to sklearn objects

"""

import contextlib
import dask.delayed
from dask.diagnostics import ProgressBar
import dask.array as da
from dask_ml.decomposition import PCA
import dask.multiprocessing
import logging
import os
import numpy as np
import plyvel
import shutil
import tables as pt
import time
import zarr
from scipy.sparse import csr_matrix
from subprocess import call

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
    kmer      = pt.StringCol(31)
    idx       = pt.UInt32Col()

class Genomes(pt.IsDescription):
    genome    = pt.StringCol(11)
    start     = pt.UInt32Col()
    stop      = pt.UInt32Col()
    idx       = pt.UInt32Col()


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

    def __init__(self, storage_dir, counter_object):

        self.counter = counter_object

        # Create storage location
        if not os.path.exists(storage_dir):
            os.makedirs(storage_dir)

        self.pytables_filepath = os.path.join(storage_dir, 'kmer_coo_pytables.h5')
        self.leveldb_filepath = os.path.join(storage_dir, 'kmer_coo_leveldb')
        self.zarr_filepath = os.path.join(storage_dir, 'kmer_dummy_zarr')

        self._tmp_kmer_counter = 0
        self._tmp_genome_counter = 0
        self.buffer = 1000000
       
    
    def enumerate(self, fasta_file_list):
        """Create pytable objects and fill with kmer data

        Overwrites any existing files

        Args:
            fasta_file_list(list): List of fasta file paths, one per line

        Returns:
            None

        """

        # Pytables
        h5file = pt.open_file(self.pytables_filepath, mode='w')
        self.init(h5file)
        if os.path.exists(self.leveldb_filepath):
            shutil.rmtree(self.leveldb_filepath)
        db = plyvel.DB(self.leveldb_filepath, create_if_missing=True)
        
        with open(fasta_file_list, 'r') as infh:
            for f in infh:
                f = f.strip()
                self.count(f, h5file, db)

        h5file.close()
        db.close()


    def init(self, h5file):
        """Create pytable objects

        Args:
            h5file(pytables.File)

        Returns:
            None

        """

        # Pytables
        h5file.create_table(h5file.root, 'indexes', COO, "COO Indexes", expectedrows=1000000000)
        h5file.create_table(h5file.root, 'kmers', Kmers, "Kmers Indexes", expectedrows=200000000)
        h5file.create_table(h5file.root, 'genomes', Genomes, "Kmers Indexes", expectedrows=2000)
        
        h5file.flush()


    def count(self, fasta_file, h5file, db):
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
            self.load(iterkmers, h5file, db, genomeid)


    def load(self, iterkmers, h5file, db, genomeid):

        wb = db.write_batch()
        kmer_memory = {}

        indextable = h5file.root.indexes
        kmertable = h5file.root.kmers
        genometable = h5file.root.genomes

        start_time = time.time()

        genome_idx = self._tmp_genome_counter
        self._tmp_genome_counter += 1
        startblock = indextable.nrows
        
        i = 0
        for kmerbstr in iterkmers:

            t0 = time.time()

            # Is kmer in memory or db
            result = None
            if kmerbstr in kmer_memory:
                result = kmer_memory[kmerbstr]

            else:
                result = db.get(kmerbstr)
                if result:
                    result = int(result)
               
            t1 = time.time()

            if not result:
                # New kmer
                kmer_idx = self._tmp_kmer_counter
                kmerrow = kmertable.row
                kmerrow['kmer'] = kmerbstr
                kmerrow['idx'] = kmer_idx
                kmerrow.append()
                self._tmp_kmer_counter += 1

                kmer_memory[kmerbstr] = kmer_idx
                wb.put(kmerbstr, str(kmer_idx).encode())

                if len(kmer_memory) >= self.buffer:
                    wb.write()
                    kmer_memory.clear()

            else:
                kmer_idx = result

            t2 = time.time()

            # Save new data point
            indexrow = indextable.row
            indexrow['row'] = genome_id
            indexrow['col'] = kmer_idx
            indexrow.append()

            t3 = time.time()

            i += 1
            if i % 100000 == 0:
                print("{} completed. Timing: query {}, insert kmer row {}, insert index row {}".format(i,
                    t1-t0, t2-t1, t3-t2))
                diff = int(t3 - start_time)
                minutes, seconds = diff // 60, diff % 60
                print('\tElapsed type ' + str(minutes) + ':' + str(seconds).zfill(2))

        endblock = indextable.nrows + i

        grow = genometable.row
        grow['start'] = startblock
        grow['stop'] = endblock
        grow['idx'] = genome_idx
        grow['genome'] = genomeid.encode()
        grow.append()

        h5file.flush()
        wb.write()

        
    def csr_matrix(self):
        """Return scipy csr_matrix representation of kmer table

        Args:
            None

        Returns tuple:
            [1]: csr_matrix
            [2]: numpy array with genome/row index
            [3]: numpy array with kmer/col index

        """

        h5file = pt.open_file(self.pytables_filepath, mode='r')

        n = h5file.root.indexes.nrows

        cm = csr_matrix(([1]*n, (h5file.root.indexes.cols.row, h5file.root.indexes.cols.col)))
        genomes = h5file.root.genomes.read(field='genome')
        kmers = h5file.root.kmers.read(field='kmer')

        h5file.close()

        return (cm, genomes, kmers)


    def dummy_matrix_sequential(self):
        """Convert COO matrix format in Pytables to zarr array with
        One-Hot-Encoding or Dummy matrix format

        BAD PERFORMANCE

        Args:
            None

        Returns:
            zarr.array

        """

        def assign(za, chunk):
            for row in chunk:
                za[row['row'],row['col']] = 1 


        h5file = pt.open_file(self.pytables_filepath, mode='r')
        ncol = h5file.root.kmers.nrows
        nrow = h5file.root.genomes.nrows
        idxtable =  h5file.root.indexes
        print(nrow, ncol)

        za = zarr.zeros((nrow, ncol), chunks=(nrow, 1000), dtype='u1')

        bsize = 1000
        t = h5file.root.indexes.nrows
        t0 = time.time()
        for b in range(0, t, bsize):
            assign(za, idxtable[b:(bsize+b)])
            if b % 100000 == 0:
                t1 = time.time()
                diff = int(t1 - t0)
                minutes, seconds = diff // 60, diff % 60
                print("{} complete.\tElapsed time {}:{}".format(round(b/t * 100, 1), str(minutes), str(seconds).zfill(2)))

        # Write to file
        # store = zarr.DirectoryStore(self.zarr_filepath)
        # group=zarr.hierarchy.group(store=store,overwrite=True,synchronizer=zarr.ThreadSynchronizer())
        # za2=group.empty('ohe',shape=za.shape,dtype=za.dtype,chunks=za.chunks)
        # za2[...]=za[...]

        print(za[:10,:10])

        return za


    def dummy_matrix(self):
        """Convert COO matrix format in Pytables to zarr array with
        One-Hot-Encoding or Dummy matrix format

        Args:
            None

        Returns:
            zarr.array

        """

        @dask.delayed
        def delayed_append(zrows):
            za = zrows[0]
            for z in zrows[1:]:
                za.append(z, axis=0)
            
            return za
           

        @dask.delayed
        def delayed_assign(chunk, ncol, j):
            zrow = zarr.zeros((1, ncol), chunks=(1,1000), dtype='u1')
            t = chunk.shape[0]
            t0 = time.time()
            i=0
            for row in chunk:
                zrow[0,row['col']] = 1
                i+=1
                if i % 100000 == 0 and j == 0:
                    t1 = time.time()
                    diff = int(t1 - t0)
                    minutes, seconds = diff // 60, diff % 60
                    print("Job {} is {} complete.\tElapsed time {}:{}".format(j,round(i/t * 100, 1), str(minutes), str(seconds).zfill(2)))

            return(zrow)

        h5file = pt.open_file(self.pytables_filepath, mode='r')
        ncol = h5file.root.kmers.nrows
        nrow = h5file.root.genomes.nrows
        idxtable =  h5file.root.indexes
        genometable =  h5file.root.genomes
        print(nrow, ncol)

        zrows = [ delayed_assign(idxtable[row['start']:row['stop']],ncol,row['start']) for row in genometable ]
        za = delayed_append(zrows)
        #za.visualize('tmp.svg')
        za = za.compute(get=dask.multiprocessing.get)
        
        # Write to file
        store = zarr.DirectoryStore(self.zarr_filepath)
        group=zarr.hierarchy.group(store=store,overwrite=True)
        za2=group.empty('dummy',shape=za.shape,dtype=za.dtype,chunks=(None, za.chunks[1]))
        za2[...]=za[...]

        h5file.close()
        
        return za


    def sample(self, p):
        """

        Args:
            None

        Returns:
            zarr.array

        """


        group = zarr.open(self.zarr_filepath, mode='r')
        darr1 = da.from_array(group['dummy'], chunks=group['dummy'].chunks)

        with ProgressBar():
            
            
            darr1 = darr1[:,darr1.sum(axis=0) > 1]
            darr1 = darr1.compute()
            ncols = darr1.shape[1]
            idx = np.random.randint(0,ncols,int(ncols*p))
            darr1 = darr1[:,idx]

            darr2 = da.from_array(darr1, chunks=(darr1.shape[0],1000))
            # darr1 = darr1.compute()
            # print(darr1)

            # # darr2 = da.from_array(darr1, chunks=darr1.chunks)
            
            # #svd_r = dd.TruncatedSVD(components=3, algorithm="tsqr", random_state=42)
            pca = PCA(n_components=3, svd_solver='randomized',random_state=0, iterated_power=4)
            # #pca = PCA(n_components=2, random_state=34, svd_solver='randomized')
            r = pca.fit(darr2)
            print(r)
            
            

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.DEBUG, format=log_fmt)


    #main(snakemake.input[0])
