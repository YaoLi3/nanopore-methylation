# ! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 28/02/2018
"""
from nanoporereads import load_sam_file, locate_snps, read_imprinted_data, get_overlapped_reads
from snps import load_VCF, map_reads, find_most_reads_snp, load_vcf_file, assemble_haplotypes
from haplotypes import *
from handlefiles import save_objects, load_objects


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
all_snps = load_objects("data/snps.obj")
reads = load_objects("data/reads.obj")
#reads_ir = load_objects("data/reads_ir.obj")

"""HMM, clustering SNPs into 2 possible haplotypes"""
# Simulation
#dummy_reads = load_objects("data/dummy/dr1.obj")  # 1000 reads  # self.bases, snps_id work
#dummy_snps = load_objects("data/dummy/ds1.obj")  # 200 snps
#dm1 = run_model(dummy_snps, dummy_reads, 10)
#dm2 = run_model(dummy_snps, dummy_reads, 10)
#compare_models(dm1, dm2)

# Real Data
m1 = run_model(all_snps, reads, 10, updateAll=False)
m2 = run_model(all_snps, reads, 10, updateAll=False)
compare_models(m1, m2)
#model = HMM(all_snps, ["Parental", "Maternal"])
#model.init_emission_matrix()
#sr_dict = model.snp_read_dict(reads)
#pm_llhd = model.assign_reads(reads)
#model.alter_update(reads, sr_dict, pm_llhd)

"""Visualize haplotypes"""
#s1, s2 = dm.read_results()
# sequence: means reference sequence
# G_Feature: one block that stands for a sequence, align to ref seq
# can add multiple G Feature blocks...
#record = GraphicRecord(sequence=s1, sequence_length=len(s1), features=[
    #GraphicFeature(start=0, end=len(s1), strand=+1, color='#ffcccc')
#])

#ax, _ = record.plot(figure_width=5)
#record.plot_sequence(ax)
#record.plot_translation(ax, (8, 23), fontdict={'weight': 'bold'})
#ax.figure.savefig('haplotype1.png', bbox_inches='tight')
