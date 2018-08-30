# ! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 08/02/2018
"""
import pysam
import time
from src.handlefiles import save_objects
from src.snps import load_VCF


#############################
#  Human Imprinted Regions  #
#############################
class ImprintedRegion:
    def __init__(self, name, chr, start, end, status=""):
        self.gene = name
        self.chr = "chr" + chr
        self.start = int(start)  # hg38
        self.end = int(end)
        self.status = status
        self.reads = []  # READS have overlapped with this imprinted region

    def __str__(self):
        return "Gene:{}, chr:{}, {}:{}".format(self.gene, self.chr, self.start, self.end)


def read_imprinted_data(fn):
    """
    Obtain human imprinted regions information.
    :param fn: string, file name
    :return: a list of ImprintedRegion objects
    """
    regions = []
    f = open(fn, "r")
    for line in f:
        line = line.strip().split("\t")
        IR = ImprintedRegion(line[0], line[3], line[1], line[2])
        regions.append(IR)
    f.close()
    return regions


############################
# Nanopore sequencing data #
############################
class NanoporeRead:
    """
    A nanopore sequencing snp instance. NA12878 data.
    """

    def __init__(self, sequencer_id, chrom, start, end, length,
                 cigar_tuples, matches_pos, quality, fastq_seq=None,
                 rs=None):
        # Basic attr
        self.id = sequencer_id
        self.chr = chrom
        self.start = start; self.end = end; self.len = length  # genomic position of first read base
        self.cigar = cigar_tuples
        self.quality = quality  # mapping quality
        self.aligned_pairs = matches_pos  # dict {ref pos: seq index}
        self.seq = fastq_seq
        self.raw_signal = rs

        # SNP attr
        self.bases = {}
        self.snps = []
        self.snps_id = []

        # Imprinted regions attr
        self.if_ir = False
        self.region = None
        self.gene_name = ""

    def if_in_imprinted_region(self, ip_regions):
        """Decide if the snp has at least one base overlapping with any imprinted region."""
        for region in ip_regions:
            if self.chr == region.chr:
                if region.start <= self.start+1 <= region.end or region.start <= self.end+1 <= region.end:
                    self.region = region
                    self.gene_name = region.gene
                    self.if_ir = True
        return self.if_ir

    def detect_snps(self, SNPs_data):
        """
        Find SNPs the read span.
        Extract allele on these SNP positions.
        """
        for snp_id, snp in enumerate(SNPs_data):
            if snp.chr == self.chr and self.start <= snp.pos-1 <= self.end:
                base = self.get_base(snp.pos-1)
                if base is not None:
                    self.bases[snp_id] = base
                    self.snps.append(snp)
                    self.snps_id.append(snp_id)

    def get_base(self, pos):
        """
        Extract sequence bases on given position.
        :param pos: ref genomic position
        :return: (str) a bases on sequence
        """
        try:
            return self.seq[self.aligned_pairs[pos]]
        except KeyError:
            pass
        except IndexError:
            pass

    def get_snp_id(self):
        """Return a list of snp_id"""
        return list(self.bases.keys())

    def get_seq(self):
        """Return sequence of the read."""
        return self.seq

    def get_base_by_snpID(self, snp_idx):
        """snp_id - base"""
        return self.bases[snp_idx]

    def set_bases(self, bases):
        """List of bases"""
        self.bases = bases

    def set_raw_signals(self, rs):
        """
        Set the raw signal for snp.
        :param rs: (numpy array) raw signal
        """
        self.raw_signal = rs

    def set_seq(self, fastq):
        """Set sequence for the read."""
        self.seq = fastq

    def __str__(self):
        """Override the default String behavior"""
        return "{}:\t{}:{}".format(self.id, self.chr, self.start)

    def __eq__(self, other):
        """Override the default Equals behavior"""
        return self.id == other.id

    def __ne__(self, other):
        """Override the default Unequal behavior"""
        return self.id != other.id or self.bases != other.bases


def load_sam_file(samfile, snps, chrom, min_qlt=10, min_len=10000):
    """
    Extract reads information from a given SAM file.
    :param samfile: (str) sam file path
    :param snps: (SNP) a list of SNP objects. SNP data on the chromosome
    :param chrom: (int/str) chromosome number
    :param min_qlt: (int) minimum mapping quality
    :param min_len: (int) minimum query length or query alignment length
    :return: a list of NanoporeRead objects
    """
    reads = []
    sf = pysam.AlignmentFile(samfile, "r")
    for read in sf:
        if read.cigartuples is not None \
                and read.reference_name == ("chr" + str(chrom)) \
                and min_qlt < read.mapping_quality \
                and min_len <= read.query_length:
            r = NanoporeRead(read.query_name,
                             read.reference_name,
                             read.reference_start,
                             read.reference_end,
                             read.infer_read_length(),
                             read.cigartuples,
                             list_to_dict(read.get_aligned_pairs(matches_only=True)),
                             read.mapping_quality,
                             fastq_seq=read.query_sequence)
            r.detect_snps(snps)
            if not r.bases == []:
                reads.append(r)
    sf.close()
    return reads


def get_overlapped_reads(reads, regions):
    """
    For a given sam file, find out if there is any READs in the file
    is located in human genetic imprinted regions.
    """
    overlapped_reads = []
    for read in reads:
        if read.if_in_imprinted_region(regions):
            overlapped_reads.append(read)
    return overlapped_reads


def locate_snps(reads, snps):
    """
    :param reads: list of Reads
    :param snps: list of SNPS
    :return list of READS that have snps in
    """
    return [read.detect_snps(snps) for read in reads]


def count_bases(reads, all_snps):
    """
    :param reads: a list of NanoporeRead objects
    :param all_snps: a list of SNP objects
    """
    ref = 0
    alt = 0
    other = 0
    for read in reads:
        for snp_id in read.bases:
            if all_snps[snp_id].ref == read.bases[snp_id]:
                ref += 1
            elif all_snps[snp_id].alt == read.bases[snp_id]:
                alt += 1
            else:
                other += 1
    print("ref+alt base num/total base num ", (ref + alt) / (ref + alt + other))
    print("other base num/total base num ", other / (ref + alt + other))


def list_to_dict(l):
    """Turn a l of tuples into a dictionary."""
    return {k: v for v, k in l}


def get_indels_regions(read):
    """Get indel region start and end positions of a read."""
    indels_blocks = []
    aligned_regions = read.get_blocks()
    start = read.reference_start; end = read.reference_end
    indels_blocks.append((start, aligned_regions[0][0]))
    for i in range(len(aligned_regions)-1):
        indels_blocks.append((aligned_regions[i][1], aligned_regions[i+1][0]))
    indels_blocks.append((aligned_regions[-1][1], end))
    return indels_blocks


if __name__ == "__main__":
    imprinted_regions = read_imprinted_data("../data/ip_gene_pos.txt")
    all_snps = load_VCF("../data/chr19.vcf")

    start = time.clock()
    reads = load_sam_file("../data/chr19_merged.sam", "19", all_snps, 10, 10000)
    elapsed = (time.clock() - start)
    print("load_sam_file: Time used:", elapsed)

    count_bases(reads, all_snps)

    save_objects("../data/snps.obj", all_snps)
    save_objects("../data/reads_no_indels.obj", reads)

    # Find READS overlapping with any human imprinted region
    o = get_overlapped_reads(reads, imprinted_regions)
    save_objects("../data/chr19_reads_ir_matched.obj", o)
