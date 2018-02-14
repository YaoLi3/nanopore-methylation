"""
__author__ = Yao LI
__email__ = liyaoo1012@163.com
__date__ = 08/02/2018
"""
import numpy as np


#######################
# Hidden Markov Model #
#######################
class HMM_haplotypes:
    """
    Clusters Nanopore reads into two haplotype groups, based on the SNPs they have (mapped to).
    This Hidden Markov Model calculates the possibility of a reads maps to parental or maternal
    reference chromosome with R(ref base) or N(alternative base) SNPs. ps. we only consider heterozygosity
    (1|0, 0|1, etc.)
    P(Ri | Hj) = -||(k) P(Ri,k | Hj,k)
    Iterate => update the possibility in the model.
    AN object is a whole model?
    """
    def __init__(self):
        self.init = [0.5, 0.5]
        self.H1 = 0
        self.H2 = 0
        self.heter = {}

    def load_vcf(self, file):
        """

        :param file:
        :return:
        """
        try:
            f = open(file, "r")
        except IOError:
            print("IOError, this vcf file is not available.")
        try:
            cnt = 0
            for line in f:
                if line.startswith("chr"):
                    chrom, pos, id, ref, alt, qual, filter, info, format, na = line.strip().split("\t")
                    h1, h2 = na.split("|")
                    if chrom == "chr19" and h1 != h2:
                        self.heter[cnt] = (chrom, pos, id, ref, alt, qual, filter, info, format, na)
                        cnt += 1
        except ValueError:
            raise RuntimeError("Not the right values to unpack.")
        return self.heter

    def cal_prob(self):
        pass


if __name__ == "__main__":
    test = HMM_haplotypes()
    tt = test.load_vcf("chr19.vcf")