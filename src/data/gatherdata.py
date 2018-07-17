"""
Gather and label data.
Should create x_train, y_train, x_test, and y_test equivalants
"""

import os
import random
import pickle
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from kmerprediction.kmer_counter import get_counts


def train_test_split(directory, percent):
    all_files = [directory + x for x in os.listdir(directory)]
    random.shuffle(all_files)
    cutoff = int(percent*len(all_files))
    train_files = all_files[:cutoff]
    test_files = all_files[cutoff:]
    return train_files, test_files


def label_data(files, metadata_sheet, genome_header, label_header):
    xl = pd.ExcelFile(metadata_sheet)
    df = xl.parse(xl.sheet_names[0])
    labels = []
    for f in files:
        label = str(df.loc[df[genome_header] == f][label_header].values[0])
        labels.append(label)
    return labels

directory = snakemake.input[0]
metadata = snakemake.input[1]
database = snakemake.input[2]
percent = float(snakemake.config["train_size"])
genome_header = snakemake.config["genome_header"]
mic_header = snakemake.config["MIC_header"]

train_files, test_files = train_test_split(directory, percent)

train_genomes = [x.replace('.fasta', '').replace(directory, '') for x in train_files]
test_genomes = [x.replace('.fasta', '').replace(directory, '') for x in test_files]

train_labels = label_data(train_genomes, metadata, genome_header, mic_header)
test_labels = label_data(test_genomes, metadata, genome_header, mic_header)

all_labels = train_labels + test_labels
le = LabelEncoder()
le.fit(all_labels)
train_labels = le.transform(train_labels)
test_labels = le.transform(test_labels)

train_data = get_counts(train_files, database)
test_data = get_counts(test_files, database)

with open(snakemake.output[0], 'wb') as f:
    pickle.dump(train_data, f)
with open(snakemake.output[1], 'wb') as f:
    pickle.dump(train_labels, f)
with open(snakemake.output[2], 'wb') as f:
    pickle.dump(test_data, f)
with open(snakemake.output[3], 'wb') as f:
    pickle.dump(test_labels, f)

with open(snakemake.output[4], 'wb') as f:
    pickle.dump(le, f)
