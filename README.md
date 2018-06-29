# nanopore-methylation
	By adapting EM algorithm in a Hidden Markov Model (HMM), we are able to cluster nanopore reads into their parental 	haplotypes. 
	One can choose use all nanopore data from the chromosome or just reads mapped to imprinted regions. 
	By analysing nanopore reads raw signals, there should be one haplotype group shows methylation happen. 
	And we are able to located methylated bases.

# Handle NA12878 data

	imprinted_regions = read_imprinted_data("data/ip_gene_pos.txt")　　
	SNPS = load_VCF("data/chr19.vcf")　　
	READS = load_sam_file("data/chr19_merged.sam", "19", SNPS) # specify chromosome number　　

Find READS overlapping with any human imprinted region

	o = get_overlapped_reads(READS, imprinted_regions)

Save pre-processed data

	save_objects("data/snps.obj", SNPS)
	save_objects("data/READS.obj", READS)
	save_objects("data/reads_ir.obj", o)

Load pre-processed data
	
	all_snps = load_objects("../data/snps.obj")
	reads = load_objects("../data/READS.obj")
	reads_ir = load_objects("data/reads_ir.obj")

# HMM clusters nanopore reads into 2 possible haplotypes
Simulation

	dummy_reads = load_objects("../data/dummy/dr1.obj")
	dummy_snps = load_objects("../data/dummy/ds1.obj")
	dm1 = run_model(dummy_snps, dummy_reads, 10)
	dm2 = run_model(dummy_snps, dummy_reads, 10)
	compare_models(dm1, dm2)

 Run model
 
	# Update all SNP sites
	model = run_model(all_snps, reads, iter_num=100) # can choose to run model iter_num times
	# Update part of SNP sites
	model = run_model(all_snps, reads, 100, updateAll=False, p=0.5) # randomly choose 50% SNP sites to update
	
	# Result: haplotypes (dict): {"m1": {snp_id: allele1}, "m2": {snp_id: allele}}
	haplotypes = model.get_haplotypes()  # each SNP site return one allele
	
	# Other results:
	# Alleles clusters saved in a dictionary {"m1": {snp_id: [allele1, allele2]}, "m2": {snp_id: [allele1, allele2]}}
	alleles_dict = model1.get_alleles()
	sorted_alleles_dict = model1.get_ordered_alleles() # sorted by SNP positions
	# haplotype strings, "-" on unassigned positions
	h1, h2 = model1.align_alleles()
	
	# Use trained model to cluster new nanopore reads data set on the same chromosome
	haplotypes_dict = model1.cluster(new_reads)
	
	# Compare model grouping results
	model2 = run_model(all_snps, reads, iter_num=100)
	compare_models(model1, model2) # return common/different alleles and SNP positions they are on
	
	# Compare results of two trained models clustering new data set
	compare_clustering(model1, model2, new_data)
