from pkg_resources import resource_filename
from sophia.utils import fastq_peek
import os

if __name__ == "__main__":
    print('Generate report')
    fastq1 = resource_filename('sophia.data', 'aln1.fastq.gz')
    fastq2 = resource_filename('sophia.data', 'aln2.fastq.gz')
    bam = resource_filename('sophia.data', 'aln.bam')

    primer_file = resource_filename('sophia.data', 'primers.txt')
    adapter_file = resource_filename('sophia.data', 'adapters.txt')
    output_dir = os.path.abspath(os.path.join(os.path.dirname(bam), '..', '..', 'output'))
    output_file = os.path.join(output_dir, 'sophia_counts.txt')

    test = False
    fastq_peek.generate_fastq_counts(fastq1, fastq2, primer_file, adapter_file, output_file, test=test)
