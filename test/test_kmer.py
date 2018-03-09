import unittest
import numpy as np


from src.data.kmer import KmerTable

class TestKmerTable(unittest.TestCase):

    def test_csrmatrix(self):

        dsk_opts = {
            '-kmer-size': 4,
            '-nb-cores': 2,
            '-abundance-min': 1,
            '-verbose': 1,
            '-out-dir': '/tmp/',
            '-out-tmp': '/tmp/'
        }

        kmer_table = KmerTable('test/tmp/dsk.h5', dsk_opts)

        kmer_table.populate('test/etc/test_genomes.txt')

        cm, genome_idx, kmer_idx = kmer_table.csr_matrix()

        expected = [b'AAAA', b'AAAT', b'AATT', b'AGGG', b'CCCC', b'CGGG', b'CCGG']
        print(kmer_idx)
        row = np.where(genome_idx == b'genome2')
        cols = np.in1d(kmer_idx, expected)
        print(cols)

        print(cm.toarray())
            

# [[1, 1, 1, 1, 1, 1, 0, 0],
#              [1, 1, 1, 1, 1, 1, 0, 0],
#              [1, 1, 1, 1, 1, 1, 1, 1]]
#         print(kmer_idx[0:])




                

if __name__ == '__main__':
    unittest.main()