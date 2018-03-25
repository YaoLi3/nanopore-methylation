# ! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 28/02/2018
"""
from nanoporereads import *
from haplotypes import *
from rawsignal import *


#############
#  Testing  #
#############
if __name__ == "__main__":
    """Nanopore reads data"""
    #DATA = NanoporeReads("data/chr19_merged.sam", "19")
    #DATA.get_reads()  # 45946 reads
    #overlap = DATA.find_imprinted(ImprintedRegions("data/ip_gene_pos.txt").get_regions(), 0, True, "data/find_imprinted_result.txt")
    #  375 reads

    """reads & snps data"""
    snps_data = load_VCF("data/chr19.vcf")  # list, 50746 SNPs on chr19
    reads_data = process_all_reads("data/find_imprinted_result.txt", snps_data)  # list, 322 reads

    """train HMM"""
    h = Hmm(snps_data, reads_data, ["A", "G", "C", "T"], ["P", "M"])
    init_random_emission = h.init_emission_matrix()
    result = minimize(h.em, init_random_emission, args=reads_data, method='BFGS')
    print(result.x)

    """raw signals"""
    #raws = get_raw_dirc("/shares/coin/yao.li/data/basecall_pass/", "/shares/coin/yao.li/signal/", overlap)
    #h1, h2 = find_haplotype(raws, haplotypes)