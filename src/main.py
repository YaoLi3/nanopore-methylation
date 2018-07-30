# ! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 28/02/2018
"""
from src.nanoporereads import *
from src.snps import *
from src.handlefiles import save_objects, load_objects

"""Load NA12878 data"""
imprinted_regions = read_imprinted_data("../data/ip_gene_pos.txt")
SNPS = load_VCF("../data/chr19.vcf")  # 50745 SNPS
READS = load_sam_file("../data/chr19_merged.sam", "19", SNPS)  # 8969 filtered READS out of 50581 total rs
# Find READS overlapping with any human imprinted region
o = get_overlapped_reads(READS, imprinted_regions)  # 86

"""Save pre-processed data"""
save_objects("../data/snps.json", SNPS)
save_objects("../data/READS.json", READS)
save_objects("../data/reads_ir.json", o)

"""Load pre-processed data"""
#all_snps = load_objects("../data/snps.obj")
#ip_snps = []
#for snp in all_snps:
    #if 10133345 <= snp.pos <= 10231286:
        #ip_snps.append(snp)
#print(ip_snps[-1].id)

#reads = load_objects("../data/READS.obj")
#reads_ir = load_objects("data/reads_ir.obj")

"""HMM, clustering SNPS into 2 possible haplotypes"""
# Simulation
#dummy_reads = load_objects("../data/dummy/dr1.obj")  # 1000 READS  # self.bases, snps_id work
#dummy_snps = load_objects("../data/dummy/ds1.obj")  # 200 snps
#dm1 = run_model(dummy_snps, dummy_reads, 10)
#print(dm1.align_alleles()[0])
#dm2 = run_model(dummy_snps, dummy_reads, 10)
#compare_models(dm1, dm2)
