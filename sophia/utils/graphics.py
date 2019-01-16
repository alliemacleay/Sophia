import numpy as np
#import matplotlib.pyplot as plt
from collections import Counter
#%matplotlib inline

def get_counts(df_comb, feature, xmax=50):
    counts = Counter(df_comb[df_comb[feature] != -1][feature].values)
    return [counts.get(i, 0) for i in range(1, xmax)]

def draw_primers_adapter(df_comb):
    xmax = df_comb.max().max() + 2
    pos = np.arange(1, xmax, 1.)
    features = {k: '' for k in ['primer_pos', 'adapter_pos', 'primer_id', 'primer_pos2', 'adapter_pos2', 'primer_id2']}
    for column in features:
        features[column] = get_counts(df_comb, column, xmax)
    primer1_ct = df_comb['primer_id'].unique
    for i in range(primers.shape[0]):
        features['primer_{}'.format(i)] = get_counts(df_comb[df_comb['primer_id'] == i], 'primer_pos', xmax)
    for i in range(primers.shape[0]):
        features['primer2_{}'.format(i)] = get_counts(df_comb[df_comb['primer_id2'] == i], 'primer_pos2', xmax)
    plt.subplot(2,1,1)
    plt.plot(pos, features['adapter_pos'], pos, features['primer_pos'])
    plt.title("Adapter and Primer positions")
    plt.ylabel("fastq 1")
    plt.subplot(2,1,2)
    plt.plot(pos, features['adapter_pos2'], pos, features['primer_pos2'])
    plt.ylabel("fastq 2")