#! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 28/02/2018
"""
from nanoporereads import *
from imprintedregions import *
from hmm import *

#############
#  Testing  #
#############
if __name__ == "__main__":
    """Nanopore reads data"""
    # DATA = NanoporeReads("data/merged.sam", "19")
    # DATA.get_reads()  # 45946 reads
    # overlap = DATA.find_imprinted(ImprintedRegions("data/ip_gene_pos.txt")
    # .get_regions(),0, True,
    # "find_imprinted_result.txt") # 375 reads

    """Human SNPs data"""
    SNP = SNPs()
    snp_data = SNP.load_VCF("data/chr19.vcf")
    snp = SNP.get_snp()  # list (pos, ref, alt, h1, h2)
    snp_reads = find_snps_in_read("find_imprinted_result.txt", snp_data)  # need to change this function

    """train HMM"""
    training, testing = split_data(snp_reads, 0.8)
    hmm = HmmHaplotypes(snp_data, training, ["Parental", "Maternal"], ["A", "T", "G", "C"])
    hmm.initialize()
    hmm.train_model(snp_reads)

    """methylation"""
