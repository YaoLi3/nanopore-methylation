# ! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 28/02/2018
"""
from nanoporereads import load_sam_file, get_snps_for_reads, read_imprinted_data, get_overlapped_reads
from snps import load_VCF, map_reads
from haplotypes import *
from handlefiles import save_objects, load_objects

#############
#  Testing  #
#############
if __name__ == "__main__":
    """Load data"""
    #imprinted_regions = read_imprinted_data("data/ip_gene_pos.txt")
    # Chr 19 data, NA12878
    #all_reads = load_sam_file("data/chr19_merged.sam", 19)  # 50581
    #all_snps = load_VCF("data/chr19.vcf")  # whole genome data

    """Link reads and SNPs"""
    #map_reads(all_snps, all_reads)  # super SLOW function
    #get_snps_for_reads(all_reads, all_snps)  # slow

    """Find reads overlapping with any human imprinted region"""
    #o = get_overlapped_reads(all_reads, imprinted_regions)  # 381

    """Save pre-processed data"""
    #save_objects("snps_data.obj", all_snps)
    #save_objects("reads_data.obj", all_reads)
    #save_objects("overlapped_reads.obj", o)

    """Load pre-processed data"""
    all_reads = load_objects("reads_data.obj")
    all_snps = load_objects("snps_data.obj")

    """HMM, clustering SNPs into 2 possible haplotypes"""
    model = Hmm(all_snps, ["P", "M"])
    iter_num = 10
    model.init_emission_matrix()

    # Whole genome analysis
    for _ in range(iter_num):
        sr_dict = model.snp_read_dict(all_reads)
        #print(sr_dict)
        pm_llhd = model.assign_reads(all_reads)
        print(pm_llhd)
        #print((np.sum(pm_llhd[:, 1]), np.sum(pm_llhd[:, 0])))
        model.update(all_reads, sr_dict, pm_llhd, hard=True, pseudo_base=1e-4)

    # Only overlapped reads analysis
