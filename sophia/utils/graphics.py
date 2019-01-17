import numpy as np
import matplotlib.pyplot as plt
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
    primer1_ct = len(df_comb['primer_id'].unique())
    primer2_ct = len(df_comb['primer_id2'].unique())
    for i in range(primer1_ct):
        features['primer_{}'.format(i)] = get_counts(df_comb[df_comb['primer_id'] == i], 'primer_pos', xmax)
    for i in range(primer2_ct):
        features['primer2_{}'.format(i)] = get_counts(df_comb[df_comb['primer_id2'] == i], 'primer_pos2', xmax)
    plt.subplot(3, 1, 1)
    plt.plot(pos, features['adapter_pos'])
    plt.title("Adapter and Primer positions")
    plt.ylabel("fastq 1 adapter")
    plt.subplot(3, 1, 2)
    plt.plot(pos, features['primer_pos'])
    plt.ylabel("fastq 1 primers")
    plt.subplot(3, 1, 3)
    plt.plot(pos, features['adapter_pos2'], pos, features['primer_pos2'])
    plt.text(5, 50, 'Adapter and primers shown together\n because adapter was not found (in blue)')
    plt.ylabel("fastq 2")
    #plt.show()
    return plt


def vcf_pie(df):
    """ create a pie chart of the variant types and the associated genes """
    def func(pct, allvals):
        """
        borrowed from the matplotlib documentation
        matplotlib.org/gallery/pie_and_polar_charts/pie_and_donut_labels.html
        """
        absolute = int(pct / 100. * np.sum(allvals))
        return "{:.1f}%\n({:d})".format(pct, absolute)

    ct = Counter(df['vartype'])
    ct2 = Counter(df[df['gene'] == 'AC073324.1']['vartype'])
    ct3 = Counter(df[df['gene'] == 'CALM1P2']['vartype'])
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(8, 6), subplot_kw=dict(aspect="equal"))
    wedges, texts, autotexts = axs[0].pie(ct.values(), autopct=lambda pct: func(pct, list(ct.values())),
                                      textprops=dict(color="w"))
    axs[0].legend(wedges, ct.keys(), title="Variant Types", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
    plt.setp(autotexts, size=8, weight="bold")
    axs[0].set_title("Variants by type")
    wedges2, texts2, autotexts2 = axs[1].pie(ct2.values(), autopct=lambda pct: func(pct, list(ct2.values())),
                                      textprops=dict(color="w"))
    axs[1].legend(wedges2, ct.keys(), title="Variant Types", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
    plt.setp(autotexts2, size=8, weight="bold")
    axs[1].set_title("AC073324.1 Variants by type")
    wedges3, texts3, autotexts3 = axs[2].pie(ct3.values(), autopct=lambda pct: func(pct, list(ct3.values())),
                                      textprops=dict(color="w"))
    axs[2].legend(wedges3, ct.keys(), title="Variant Types", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
    plt.setp(autotexts3, size=8, weight="bold")
    axs[2].set_title("CALM1P2 Variants by type")
    #plt.show()
    return plt

