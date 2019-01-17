import pandas as pd
import vcf

class Vcf(object):
    """ Class for a VCF file """
    def __init__(self, vcf_path):
        self.path = vcf_path
        self.header_rows = 0
        self.df = self._vcf_df()

    def _vcf_df(self):
        vcf_dict = {'chrom': [], 'pos': [], 'ref': [], 'alt': [], 'dp': [], 'af': [], 'vartype': [], 'mod': [], 'gene': []}
        with open(self.path, 'r') as fh:
            for rec in vcf.Reader(fh):
                vcf_dict['chrom'].append(rec.CHROM)
                vcf_dict['pos'].append(rec.POS)
                vcf_dict['ref'].append(rec.REF)
                vcf_dict['alt'].append(rec.ALT[0])

                vcf_dict['dp'].append(rec.INFO['DP'][0] + rec.INFO['DP'][1])
                vcf_dict['af'].append(rec.INFO['AF'][0])
                csq = rec.INFO['CSQ'][0].split('|')
                vcf_dict['vartype'].append(csq[1])
                vcf_dict['mod'].append(csq[2])
                vcf_dict['gene'].append(csq[3])

        df = pd.DataFrame(vcf_dict)
        return df

    def _set_header_rows(self):
        with open(self.path) as fh:
            for line in fh:
                if line.startswith("#"):
                    self.header_rows += 1
                else:
                    self.header_rows -= 1
                    break

