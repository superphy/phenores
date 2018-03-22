import unittest
import numpy as np
import zarr


from src.data.kmer import KmerTable
from src.data.kmer import DskCounter

class TestKmerTable(unittest.TestCase):

    def test_csrmatrix(self):

        dsk_opts = {
            '-kmer-size': 5,
            '-nb-cores': 2,
            '-abundance-min': 1,
            '-verbose': 1,
            '-out-dir': '/tmp/',
            '-out-tmp': '/tmp/'
        }

        dc = DskCounter(dsk_opts)
        kmer_table = KmerTable('test/tmp/', dc)
        kmer_table.buffer = 2

        kmer_table.enumerate('test/etc/test_genomes.txt')

        cm, genome_idx, kmer_idx = kmer_table.csr_matrix()

        cmarray = cm.toarray()

        expected = [b'AAAAA', b'CCCCC']
        row = np.where(genome_idx == b'genome3')
        cols = np.in1d(kmer_idx, expected)

        self.assertTrue(np.sum(cols) == 2)
        self.assertTrue(np.sum(cmarray[row,:]) == 2)
        self.assertTrue(np.sum(cmarray[row,cols]) == 2)


    def test_dummy_matrix(self):

        dsk_opts = {
            '-kmer-size': 5,
            '-nb-cores': 2,
            '-abundance-min': 1,
            '-verbose': 1,
            '-out-dir': '/tmp/',
            '-out-tmp': '/tmp/'
        }

        dc = DskCounter(dsk_opts)
        kmer_table = KmerTable('test/tmp/', dc)
        kmer_table.buffer = 2

        kmer_table.enumerate('test/etc/test_genomes.txt')
        dm = kmer_table.dummy_matrix()


        group = zarr.open('test/tmp/kmer_dummy_zarr', mode='r')
        dm2 = group['dummy']

        
        print(dm[:])
        print(dm.shape)
        print(dm.chunks)
        print(dm2[:])
        print(dm2.shape)
        print(dm2.chunks)
        self.assertIsNotNone(dm)

        #print(dm[:])
        
        # expected = [b'AAAAA', b'CCCCC']
        # row = np.where(genome_idx == b'genome3')
        # cols = np.in1d(kmer_idx, expected)

        # self.assertTrue(np.sum(cols) == 2)
        # self.assertTrue(np.sum(cmarray[row,:]) == 2)
        # self.assertTrue(np.sum(cmarray[row,cols]) == 2)
        
                

if __name__ == '__main__':
    unittest.main()