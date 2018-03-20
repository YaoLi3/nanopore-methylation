# ! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 28/02/2018
"""
from nanoporereads import *
from imprintedregions import *
from hmm import *
from rawsignal import *
from handlefiles import *

#############
#  Testing  #
#############
if __name__ == "__main__":
    """Nanopore reads data"""
    #DATA = NanoporeReads("data/chr19_merged.sam", "19")
    #DATA.get_reads()  # 45946 reads
    #overlap = DATA.find_imprinted(ImprintedRegions("data/ip_gene_pos.txt").get_regions(), 0, True, "find_imprinted_result.txt")
    #  375 reads

    """read snps data"""
    snp_data = load_VCF("data/chr19.vcf") # list, 50746 SNPs on chr19
    read_snp = process_all_reads("data/find_imprinted_result.txt", snp_data) # list, 322 reads

    """train HMM"""
    hmm = HmmHaplotypes(snp_data, read_snp, ["P", "M"], ["A", "T", "G", "C"])
    print(hmm.old_emission[0, 0,])
    print(hmm.emission[0, 0,])
    hmm.initialize()
    print(hmm.old_emission[0, 0, ])
    print(hmm.emission[0, 0, ])
    hmm.assign_reads()
    print(len(hmm.old_d0))
    print(len(hmm.d0))
    o, e = hmm.cal_emission()
    snp_index = snp_data.index(hmm.d0[0].snps[0])
    print(o[0, 0, ])
    print(e[0, 0, ])
    print(o[snp_index, 0, ])
    print(e[snp_index, 0, ])
    hmm.assign_reads()
    print(len(hmm.old_d0))
    print(len(hmm.d0))
    hmm.cal_emission()
    hmm.assign_reads()
    snp_index = snp_data.index(hmm.d0[0].snps[0])
    print(o[snp_index, 0,])
    print(e[snp_index, 0,])
    print(len(hmm.old_d0))
    print(len(hmm.d0))

    """raw signal data"""
    #raw_001 = get_raw_dirc("/shares/coin/yao.li/data/basecall_pass/", "/shares/coin/yao.li/signal/", overlap)
    #raw_000 = get_raw_dirc("/shares/coin/yao.li/data/basecall_pass/", "/shares/coin/yao.li/signal/", overlap, basecall_group="Basecall_1D_000")
    #print(np.load("data/chr19_merge.npy")[0])

    """methylation clustering"""