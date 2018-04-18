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

    """Chr 19"""
    all_reads = load_sam_file("data/chr19_merged.sam", 19)  # 50581
    o = get_overlapped_reads(all_reads, imprinted_regions)  # 381

    all_snps = load_VCF("data/chr19.vcf")  # whole genome data
    map_reads(all_snps, all_reads)  # super SLOW function TODO: optimize

    model = Hmm(all_snps, ["P", "M"])
    iter_num = 10
    model.init_emission_matrix()
    # Whole genome analysis
    for _ in range(iter_num):
        sr_dict = model.snp_read_dict(all_reads)
        pm_llhd = model.assign_reads(all_reads)
        print((np.sum(pm_llhd[:, 1]), np.sum(pm_llhd[:, 0])))
        model.update(all_reads, sr_dict, pm_llhd, hard=False, pseudo_base=1e-4)