from pkg_resources import resource_filename
from sophia.utils import fastq_peek, graphics, vcf_peek
import os
import pandas as pd

if __name__ == "__main__":
    print('Generate report')
    fastq1 = resource_filename('sophia.data', 'aln1.fastq.gz')
    fastq2 = resource_filename('sophia.data', 'aln2.fastq.gz')
    bam = resource_filename('sophia.data', 'aln.bam')
    vcf_path = resource_filename('sophia.data', 'aln.vep.annotated.vcf')

    primer_file = resource_filename('sophia.data', 'primers.txt')
    adapter_file = resource_filename('sophia.data', 'adapters.txt')

    output_dir = os.path.abspath(os.path.join(os.path.dirname(bam), '..', '..', 'output'))
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    output_fq_file = os.path.join(output_dir, 'sophia_fastq_counts.txt')
    output_fq_png = os.path.join(output_dir, 'sophia_fastq_graph.png')
    output_vcf_file = os.path.join(output_dir, 'sohpia_vcf_counts.txt')
    output_vcf_png = os.path.join(output_dir, 'sohpia_vcf_pie.png')

    test = False

    df = fastq_peek.generate_fastq_counts(fastq1, fastq2, primer_file, adapter_file, test=test)
    df.to_csv(output_fq_file, sep='\t', index=False)
    print('FASTQ output saved at {}'.format(output_fq_file))

    # df = pd.read_csv(output_fq_file, sep='\t')
    plt1 = graphics.draw_primers_adapter(df)
    # plt1.show()
    plt1.savefig(output_fq_png)
    print('FASTQ plot saved at {}'.format(output_fq_png))

    vcf = vcf_peek.Vcf(vcf_path)
    vcf.df.to_csv(output_vcf_file, sep='\t', index=False)
    print('VCF output saved at {}'.format(output_vcf_file))

    plt2 = graphics.vcf_pie(vcf.df)
    plt2.savefig(output_vcf_png)
    print('VCF plot saved at {}'.format(output_vcf_png))



