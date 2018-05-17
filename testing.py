# ! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 28/02/2018
"""
from nanoporereads import load_sam_file, locate_snps, read_imprinted_data, get_overlapped_reads
from snps import load_VCF, map_reads, find_most_reads_snp, load_vcf_file
from haplotypes import *
from handlefiles import save_objects, load_objects

#############
#  Testing  #
#############
if __name__ == "__main__":
    """Load NA12878 data"""
    #imprinted_regions = read_imprinted_data("data/ip_gene_pos.txt")
    #all_snps = load_VCF("data/chr19.vcf")  # 50745 SNPs
    #reads = load_sam_file("data/chr19_merged.sam", "19", all_snps)  # 8969 filtered reads out of 50581 total rs
    # Find reads overlapping with any human imprinted region
    #o = get_overlapped_reads(reads, imprinted_regions)  # 86

    """Save pre-processed data"""
    #save_objects("data/snps.obj", all_snps)
    #save_objects("data/reads.obj", reads)
    #save_objects("data/reads_ir.obj", o)

    """Load pre-processed data"""
    #all_snps = load_objects("data/snps.obj")
    #reads = load_objects("data/reads.obj")
    #reads_ir = load_objects("data/reads_ir.obj")

    """HMM, clustering SNPs into 2 possible haplotypes"""
    # Simulation
    dummy_reads = load_objects("dummy_reads.obj")
    dummy_snps = load_objects("dummy_snps.obj")
    #map_reads(all_snps, reads)
    #for read in reads:
        #for snp in read.snps:
            #print(snp)

    m1 = run_model(dummy_snps, dummy_reads, 10, False)
    #m2 = run_model(dummy_snps, dummy_reads, 10, False)
    #compare_assignments(m1.get_snps(), m2.get_snps())
    #gold_standard(dummy_snps, m1.get_snps())
