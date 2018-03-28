# ! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 28/02/2018
"""
from nanoporereads import *
from haplotypes import *
from rawsignal import *
from snps import *

#############
#  Testing  #
#############
if __name__ == "__main__":
    """Nanopore reads data"""
    # extract_fastq("/dict/fast5_folder/")
    # Use Minimap2 map reads to reference, get SAM file
    # DATA = NanoporeReads("data/chr19_merged.sam", "19")
    # DATA.get_reads()  # 45946 reads
    # overlap = DATA.find_imprinted(ImprintedRegions("data/ip_gene_pos.txt").get_regions(), 0, True, "data/find_imprinted_result.txt")
    #  375 reads

    """reads & snps data"""
    snps_data = load_VCF("data/chr19.vcf")  # list, 50746 SNPs on chr19
    reads_data = process_all_reads("data/find_imprinted_result.txt", snps_data)  # list, 302 reads

    """train HMM"""
    model = Hmm(snps_data, reads_data, ["A", "G", "C", "T"], ["P", "M"])
    # initial
    matrix_init = model.init_emission_matrix()
    m0, m1, m0_pos, m1_pos = model.assign_reads(matrix_init)
    # 1st round
    model.maximize_likelihood_for_each_snp_pos(matrix_init, 0, m0, m0_pos)  # NOTE: m0_pos right now is np array, zeros, float64. does not have chrom attribute
    matrix_old = model.maximize_likelihood_for_each_snp_pos(matrix_init, 1, m0, m0_pos)
    m01, m11, m0_pos1, m1_pos1 = model.assign_reads(matrix_old)
    # 2nd round
    model.maximize_likelihood_for_each_snp_pos(matrix_old, 0, m01, m0_pos1)
    matrix_new = model.maximize_likelihood_for_each_snp_pos(matrix_old, 1, m01, m0_pos1)
    m02, m12, m0_pos2, m1_pos2 = model.assign_reads(matrix_new)

    print(len(m0), len(m1))
    print(len(m01), len(m11))
    print(len(m02), len(m12))

    """raw signals"""
    # raws = get_raw_dirc("/shares/coin/yao.li/data/basecall_pass/", "/shares/coin/yao.li/signal/", overlap)
    # h1, h2 = find_haplotype(raws, m02, m12)