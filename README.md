# nanopore-methylation

# Handle NA12878 data

	imprinted_regions = read_imprinted_data("data/ip_gene_pos.txt")　　
	SNPS = load_VCF("data/chr19.vcf")　　
	READS = load_sam_file("data/chr19_merged.sam", "19", SNPS)　　

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

HMM clustering SNPS into 2 possible haplotypes

# Simulation
	dummy_reads = load_objects("../data/dummy/dr1.obj")
	dummy_snps = load_objects("../data/dummy/ds1.obj")
	dm1 = run_model(dummy_snps, dummy_reads, 10)
	print(dm1.read_results()[0])
	dm2 = run_model(dummy_snps, dummy_reads, 10)
	compare_models(dm1, dm2)

  # Run model
  
	model = run_model(all_snps, reads, iter_num=100) # can choose to run model iter_num times
	haplotype1, haplotype2 = model.read_results()
	m2 = run_model(all_snps, reads, iter_num=100)
	compare_models(model, m1)
