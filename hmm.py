"""
__author__ = Yao LI
__email__ = liyaoo1012@163.com
__date__ = 08/02/2018
"""
import numpy as np


#######################
# Hidden Markov Model #
#######################
class HmmHaplotypes:
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

    def cal_prob(self):
        pass


def loadVCF(file):
    """
    Read a VCF file and return he data in it.
    :param file:
    :return: heter: (dictionary) data of the variants (SNPs), key are numbers
    """
    try:
        f = open(file, "r")
        heter = {}
        cnt = 0
        for line in f:
            if line.startswith("chr"):
                chrom, pos, id, ref, alt, qual, filter, info, format, na = line.strip().split("\t")
                h1, h2 = na.split("|")
                if chrom == "chr19" and h1 != h2:
                    heter[cnt] = (chrom, pos, id, ref, alt, qual, filter, info, format, h1, h2)
                    cnt += 1
        f.close()
        return heter
    except ValueError:
        raise RuntimeError("Not the right values to unpack.")
    except IOError:
        raise IOError("This vcf file is not available.")


if __name__ == "__main__":
    snp_data = loadVCF("data/chr19.vcf")
    for i in snp_data:
        print(snp_data[i])
