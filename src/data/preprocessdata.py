from kmerprediction.feature_scaling import scale_to_range
from kmerprediction.feature_selection import select_k_best
import pickle 

with open(snakemake.input[0], 'rb') as f:
    train_data = pickle.load(f)
with open(snakemake.input[1], 'rb') as f:
    train_labels = pickle.load(f)
with open(snakemake.input[2], 'rb') as f:
    test_data = pickle.load(f)
with open(snakemake.input[3], 'rb') as f:
    test_labels = pickle.load(f)

data = [train_data, train_labels, test_data, test_labels]
data, _ = select_k_best(data, None, k=270)
data = scale_to_range(data)

with open(snakemake.output[0], 'wb') as f:
    pickle.dump(data[0], f)
with open(snakemake.output[1], 'wb') as f:
    pickle.dump(data[1], f)
with open(snakemake.output[2], 'wb') as f:
    pickle.dump(data[2], f)
with open(snakemake.output[3], 'wb') as f:
    pickle.dump(data[3], f)


