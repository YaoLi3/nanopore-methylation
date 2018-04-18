# ! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 28/02/2018
"""
from nanoporereads import *
from snps import *
from haplotypes import *
from rawsignal import *

#############
#  Testing  #
#############
if __name__ == "__main__":
    imprinted_regions = read_imprinted_data("data/ip_gene_pos.txt")
    all_reads = load_sam_file("data/chr19_merged.sam", 19)  # 50581
    oreads = get_overlapped_reads(all_reads, imprinted_regions)  # 381

    all_snps = load_VCF("data/chr19.vcf", all_reads)  # whole genome data
    for snp in all_snps:
        print(snp.reads)