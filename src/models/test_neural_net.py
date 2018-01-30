import pickle
from keras.utils import to_categorical
from keras.models import load_model

model = load_model(snakemake.input[0])

with open(snakemake.input[1], 'rb') as f:
    test_data = pickle.load(f)
with open(snakemake.input[2], 'rb') as f:
    test_labels = pickle.load(f)

test_labels = to_categorical(test_labels)

score = model.evaluate(test_data, test_labels)

output = dict(zip(model.metrics_names, score))

print(output)
