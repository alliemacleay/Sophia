#!pip install pysam
# version 0.15.1
# !pip install biopython  # version 1.73

from Bio import SeqIO
import pandas as pd
import gzip


# define funtion to return a pandas dataframe with all fastq read sequences
def get_fastq_df(path):
    dct = {'seq': [], 'id': [], 'qual': []}
    with gzip.open(path, "rt") as fastq_handle:
        for record in SeqIO.parse(fastq_handle, "fastq"):
            dct['seq'].append(str(record.seq))
            dct['id'].append(record.id)
            dct['qual'].append("")
    return pd.DataFrame(dct)


# look for exact match to sequence
def find_exact(match, target):
    return target.find(match)


# find fucntion for dataframe row.  mstring can be a sequence or a list of sequences
# returns the position of the read found
def df_which_primer(row, prlist, method=find_exact):
    seq = row['seq']
    for num, single in enumerate(prlist):
        pos = method(single, seq)
        if pos != -1:
            return num
    return -1


# find fucntion for dataframe row.  prlist can be a sequence or a list of sequences
# returns the id of the primer found
def df_seq_pos(row, prlist, method=find_exact):
    seq = row['seq']
    if isinstance(prlist, list):
        for num, single in enumerate(prlist):
            pos = method(single, seq)
            if pos != -1:
                return pos
        return -1
    else:
        return method(prlist, seq)

def generate_fastq_counts(fastq1, fastq2, primer_file, adapter_file, output_file, test=False):
    # load primers and adapters
    primers = pd.read_csv(primer_file, sep="\t")
    adapters = pd.read_csv(adapter_file, sep="\t", header=None, index_col=0)
    r1_primers = list(primers["forward"].values)
    r1_adapter = adapters.loc["Adaptor1"][1]
    r2_primers = list(primers["reverse"].values)
    r2_adapter = adapters.loc["Adaptor2"][1]

    # Get fastq dataframes
    df1 = get_fastq_df(fastq1)
    df2 = get_fastq_df(fastq2)

    # subset for testing
    if test:
        df1 = df1.sample(5000)
        mask = df1.index
        df2 = df2.loc[mask]

    # find primers
    df1['primer_id'] = df1.apply(df_which_primer, args=(r1_primers,), axis=1)
    df1['primer_pos'] = df1.apply(df_seq_pos, args=(r1_primers,), axis=1)
    df1['adapter_pos'] = df1.apply(df_seq_pos, args=(r1_adapter,), axis=1)

    df2['primer_id'] = df2.apply(df_which_primer, args=(r2_primers,), axis=1)
    df2['primer_pos'] = df2.apply(df_seq_pos, args=(r2_primers,), axis=1)
    df2['adapter_pos'] = df2.apply(df_seq_pos, args=(r2_adapter,), axis=1)

    df2 = df2.rename(columns={'primer_pos': 'primer_pos2', 'adapter_pos': 'adapter_pos2', 'primer_id': 'primer_id2'})
    df_comb = pd.concat(
        [df1[['primer_pos', 'adapter_pos', 'primer_id']], df2[['primer_pos2', 'adapter_pos2', 'primer_id2']]], axis=1)
    df_comb.write_csv(output_file, sep='\t', index=False)
