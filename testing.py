# ! /usr/bin/env python
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
    #DATA = NanoporeReads("data/merged.sam", "19")
    #DATA.get_reads()  # 45946 reads
    #overlap = DATA.find_imprinted(ImprintedRegions("data/ip_gene_pos.txt").get_regions(), 0, True, "find_imprinted_result.txt")
    #  375 reads

    """read snps data"""
    snp_data = load_VCF("data/chr19.vcf") # list, 50746 SNPs
    read_snp = process_all_reads("find_imprinted_result.txt", snp_data) # list, 322 reads

    """train HMM"""
    #training, testing = split_data(read_snp, 0.8)
    hmm = HmmHaplotypes(snp_data, read_snp, ["P", "M"], ["A", "T", "G", "C"])
    #print(hmm.emission.shape) # (50746, 2, 4)
    #hmm.init_emission()
    #print(hmm.emission.shape) #(50746, 2, 4)
    hmm.initialize()
    #print(hmm.transition.shape) # (2, 2)
    #print(len(hmm.M)) # list of 161 OverlapRead object
    #hmm.cal_read_prob(read_snp[0], 0)
    hmm.cal_n_assign(read_snp)
    print(len(hmm.d0)) #131 M reads 149 155
    print(len(hmm.d1)) #163 P reads 153 147
    print(hmm.emission.shape) #(50746, 2, 4)

    """raw signal data"""

    """methylation clustering"""
