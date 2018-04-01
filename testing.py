# ! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 28/02/2018
"""
#from nanoporereads import *
#from haplotypes import *
#from rawsignal import *
#from snps import *
import math
import numpy as np
from scipy.optimize import minimize
from snps import load_VCF
from snps import process_all_reads
from haplotypes import Hmm

#############
#  Testing  #
#############
if __name__ == "__main__":
    """Nanopore reads data"""
    #extract_fastq("/dict/fast5_folder/")
    # Use Minimap2 map reads to reference, get SAM file
    # DATA = NanoporeReads("data/chr19_merged.sam", "19")
    # DATA.get_reads()  # 45946 reads
    # overlap = DATA.find_imprinted(ImprintedRegions("data/ip_gene_pos.txt").get_regions(), 0, True, "data/find_imprinted_result.txt")
    #  375 reads

    """reads & snps data"""
    snps_data = load_VCF("data/chr19.vcf")  # list, 50746 SNPs on chr19
    reads_data = process_all_reads("data/find_imprinted_result.txt", snps_data)  # list, 302 reads

    # snps = []
    # for read in reads_data:
    # for snp in read.snps:
    # if snp not in snps:
    # snps.append(snp)

    """train HMM"""


    h = Hmm(snps_data, ["P", "M"])
    iter_num = 10
    h.init_emission_matrix()
    for _ in range(iter_num):
        sr_dict = h.snp_read_dict(reads_data)
        pm_llhd = h.assign_reads(reads_data)
#        print(h.theta[9519,:])
        print((np.sum(pm_llhd[:,1]),np.sum(pm_llhd[:,0])))
        h.update(reads_data,sr_dict, pm_llhd,hard = False,pseudo_base = 1e-4)
    
    # try:
#    matrix = h.init_emission_matrix()
#    m0, m1, m0_pos, m1_pos = assign_reads(matrix, reads_data, snps_data)
#    print(len(m0),len(m1))
#    for i in range(iter_num):
#        matrix = maximize_likelihood_for_each_snp(matrix,0,m0,m0_pos,snps_data)
#        matrix = maximize_likelihood_for_each_snp(matrix,1,m1,m1_pos,snps_data)
#        m0,m1,m0_pos,m1_pos = assign_reads(matrix,reads_data,snps_data)
#        print(len(m0),len(m1))





