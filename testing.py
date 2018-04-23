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
    #all_snps = load_VCF("data/chr19.vcf")  # 50745

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
    overlap_reads = load_objects("overlapped_reads.obj")

    """HMM, clustering SNPs into 2 possible haplotypes"""
    iter_num = 100
    total_snps_assignments = np.zeros((len(all_snps), 2))
    for _ in range(iter_num):
        model_snps_assign = model_iterations(all_snps, all_reads, 10)
        if total_snps_assignments.shape == model_snps_assign.shape:
            total_snps_assignments += model_snps_assign
        else:
            raise ValueError("operands could not be broadcast together with different shapes.")

    with open("model_result.txt", "w") as f:
        f.write("Chr19 whole genome data, Oxford Nanopore reads. Iterate 100 times.\n")
        f.write("CHR\t\tPOS\t\tREF\t\tALT\t\tHMM\n")
        for snp_id in range(total_snps_assignments.shape[0]):
            m1p = total_snps_assignments[snp_id, 0] / 1000
            m2p = total_snps_assignments[snp_id, 1] / 1000
            #f.writelines("SNP assignments: model1:{}\tmodel2:{}\iter_num".format(m1p, m2p))
            snp = all_snps[snp_id]
            f.write("{}\t\t{}\t\t{}\t\t{}\t\t{}/{}\n".format(snp.chrom, snp.pos, snp.ref, snp.alt, m1p, m2p))
