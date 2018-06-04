
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
        #"data/interim/streptomycin/kmers/kmer_matrix.npz",
        "data/interim/kmers/kmer_matrix.npz",


# rule streptokmers:
#     input:
#         "data/interim/streptomycin_fasta_files.txt"
#     output:
#         dir="data/interim/streptomycin/kmers/"
#     params:
#         k=11,
#         mins=1,
#         v=0,
#         cores=12,
#         d="data/interim/dsk"
#     script: 
#         "src/data/make_kmer_table.py"


rule kmers:
    input:
        "data/interim/fasta_files_list.txt"
    output:
        "data/interim/kmers/kmer_ohe_zarr",
        "data/interim/kmers/kmer_index_leveldb"
    params:
        k=11,
        mins=1,
        v=0,
        cores=12,
        d="data/interim/dsk"
    script: 
        "src/data/make_kmer_table.py"


# rule uniquify:
#     input:
#         "data/interim/streptomycin/kmers/kmer_ohe_zarr",
#         "data/interim/streptomycin/kmers/kmer_index_leveldb"
#     output:
#         "data/interim/streptomycin/kmers/kmer_patterns_leveldb",
#         "data/interim/streptomycin/kmers/kmer_matrix.npz"
#     run:
#         import zarr
#         import pandas as pd
#         import dask
#         import dask.array as da
#         import dask.dataframe as dd
#         import mmh3
#         import numpy as np
#         import plyvel
#         import sys
#         import time

#         start_time = time.time()

#         store = zarr.DirectoryStore(input[0])
#         group = zarr.hierarchy.group(store=store,overwrite=False)
#         za = group['ohe']

#         db = plyvel.DB(output[1], create_if_missing=True)

#         nkmers = za.shape[1]
#         ngenomes = za.shape[0]

#         class Counter:
#             c = 0

#         def hashfunc(x, ngenomes):
#             arrstr = np.array2string(x, threshold=ngenomes)
#             hashstr = str(mmh3.hash(arrstr))
#             return hashstr

#         def loadfunc(x):
#             istr = str(Counter.c).encode()
#             xstr = str(x).encode()
#             if not db.get(xstr):
#                 db.put(xstr, istr)
#             db.put(istr, xstr)

#             Counter.c += 1
            
#         loadvecfunc = np.vectorize(loadfunc)      

#         stepsize = 1000
#         for i in range(0, 3000, stepsize):
#             t0 = time.time()
#             #hashstr = mmh3.hash(np.array2string(a[:,i], threshold=ngenomes))
#             selection = np.arange(i,i+stepsize)
#             na = za.get_orthogonal_selection((slice(None), selection))
#             hashstrs = np.apply_along_axis(hashfunc, 0, na, ngenomes)
#             t1 = time.time()
#             loadvecfunc(hashstrs)
#             t2 = time.time()

#             if i % stepsize == 0:
#                 print("{} completed. Timing: hashing {}, loading {}".format(i,
#                     t1-t0, t2-t1))
#                 diff = int(t1 - start_time)
#                 minutes, seconds = diff // 60, diff % 60
#                 print('\tElapsed type ' + str(minutes) + ':' + str(seconds).zfill(2))
#                 # print(hashstr.shape)
#                 # print(hashstr[0])
#                 # print(hashstr[1])


#         print(sys.getsizeof(hashstr))



rule filter:
    input:
        "data/interim/kmers/kmer_ohe_zarr",
        "data/interim/fasta_files_list.txt",
        "data/interim/kmers/kmer_index_leveldb"
    params:
        maj=2526,
        min=25,
        k=11
    output:
        "data/interim/kmers/kmer_matrix.npz"
    run:
        import zarr
        import pandas as pd
        import dask
        import dask.array as da
        import dask.dataframe as dd
        import numpy as np
        import plyvel
        from os.path import splitext, basename

        # Output genome order
        with open(input[1], 'r') as infh:
            fasta_files = infh.read().splitlines()

        genomes = [ splitext(basename(f))[0] for f in fasta_files]
        
        garray = np.array(genomes)

        print("Genome array built")

        store = zarr.DirectoryStore(input[0])
        group = zarr.hierarchy.group(store=store,overwrite=False)
        za = group['ohe']

        a = da.from_array(za, chunks=za.chunks)
        nkmers = a.shape[1]

        # Record kmer order
        db = plyvel.DB(input[2])
        dt = 'U'+str(params['k'])
        ko = np.empty(nkmers, dtype=dt)

        for k, v in db:
            # print(k)
            # print(v)
            i = int(v.decode('utf-8'))
            ko[i] = k.decode('utf-8')

        print("Kmer order array built")

        print("Start of filtering")
        with dask.set_options(scheduler='threads'):

            # Filter by counts & unique patterns
            frequencies = a.sum(axis=0)
            mask = (frequencies > params['min']) & (frequencies < params['maj'])
            a = a[:,mask]
            karray = a.compute()
            print(karray.shape)

            koarray = ko[mask]

            # sample = np.random.choice(a.shape[0],size=int(a.shape[0]*.01))
            # a2 = a[sample,:]

            # mask = da.apply_along_axis(pass_cond, 0, a).compute()
            # print(mask.shape)
            # print(np.sum(mask))
            
            # mask.compute()
            # a = a[:,mask]
            

        print("End of filtering")
        np.savez(output[0], kmers=karray, kmer_order=koarray, genome_order=garray)



            


        













