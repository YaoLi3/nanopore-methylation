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
    DATA = NanoporeReads("data/merged.sam", "19")
    DATA.get_reads()  # 45946 reads
    o = DATA.find_imprinted(ImprintedRegions
                            ("data/mart_export.txt").get_regions(),
                            0, False, "find_imprinted_result.txt")

    snp = process_data(load_VCF("data/chr19.vcf"))
    h1, h2 = split_data(split_data(snp, 0.8)[0], 0.5)
    data = find_read_snp("find_imprinted_result.txt", snp)