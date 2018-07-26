

#################################################################################
# FUNCTIONS                                                                     #
#################################################################################


#################################################################################
# GLOBALS                                                                       #
#################################################################################


#################################################################################
# RULES                                                                         #
#################################################################################


rule all:
    input:
        "data/interim/kmers/kmer_matrix.npz"


rule kmers:
    input:
        "data/interim/fasta_files_list_no_ecoli.txt"
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


rule filter:
    input:
        "data/interim/kmers/kmer_ohe_zarr",
        "data/interim/fasta_files_list_no_ecoli.txt",
        "data/interim/kmers/kmer_index_leveldb"
    params:
        maj=2260,
        min=0,
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

            # Filter by counts
            frequencies = (a != 0).sum(0)
            mask = (frequencies > params['min']) & (frequencies < params['maj'])
            a = a[:,mask]
            karray = a.compute()
            print(karray.shape)

            koarray = ko[mask]
            

        print("End of filtering")
        np.savez(output[0], kmers=karray, kmer_order=koarray, genome_order=garray)



            


        













