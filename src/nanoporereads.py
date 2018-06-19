# ! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 08/02/2018
"""
import pysam
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

    def __init__(self, sequencer_id, chrom, start, end, quality, fastq_seq=None, rs=None):
        # Basic attr
        self.id = sequencer_id
        self.chr = chrom
        self.start = int(start)
        self.end = int(end)
        self.quality = quality  # mapping quality

        self.seq = fastq_seq
        self.raw_signal = rs

        # SNP attr
        self.snps = []
        self.snps_id = []
        self.bases = {}  # bases of read on SNPS loci. a mini haplotype  # bases, not base. {snp_id: read base}
        self.gt = ""

        # Imprinted regions attr
        self.if_ir = False
        self.region = None
        self.gene_name = ""
        if self.seq is not None:
            self.segment_pos = [0, len(self.seq)]  # overlap segment, snp_id on snp
        self.seg_ref_pos = [0, 0]  # segment pos on ref genome

    def if_in_imprinted_region(self, ip_regions):
        """Decide if the snp has at least one base overlapping with any imprinted region."""
        for region in ip_regions:
            if self.chr == region.chr:
                if region.start <= self.start <= region.end or region.start <= self.end <= region.end:
                    self.region = region
                    self.gene_name = region.gene
                    self.seg_ref_pos = [self.start, self.end]
                    self.if_ir = True
        return self.if_ir

    def detect_snps(self, SNPs_data):
        """
        Find snps in one READs.
        :param SNPs_data: (list) all SNPS data
        """
        for snp_id, snp in enumerate(SNPs_data):
            if snp.chr == self.chr and self.start <= snp.pos < self.end:  # check
                self.snps.append(snp)
                self.snps_id.append(snp_id)  # not necessary
                if self.seq is not None:
                    self.bases[snp_id] = self.get_base(snp.pos)

    def get_base(self, pos):
        """
        Extract sequence bases on given position.
        :param pos: ref genome position
        :return: (str) a bases on sequence
        """
        try:
            index = int(pos) - int(self.start) - 1
            return self.seq[index]
        except IndexError:
            print(self.start, self.end, pos)

    def get_dict_base(self, snp_idx):
        return self.bases[snp_idx]

    def set_bases(self, bases):
        """List of bases"""
        self.bases = bases

    def set_raw_signal(self, rs):
        """
        Set the raw signal for snp.
        :param rs: (numpy array) raw signal
        """
        self.raw_signal = rs

    def set_seq(self, fastq):
        """Set sequence for the read."""
        self.seq = fastq

    def __str__(self):
        return "{}:\t{}:{},{}".format(self.id, self.chr, self.start, self.end)

    def __eq__(self, other):
        """Override the default Equals behavior"""
        return self.id == other.id

    def __ne__(self, other):
        """Override the default Unequal behavior"""
        return self.id != other.id or self.snps != other.snps or self.start != other.start


def load_sam_file(samfile, chr, snps):
    reads = []
    sf = pysam.AlignmentFile(samfile, "r")
    for read in sf:
        if read.reference_name == ("chr" + chr) and 10 < read.mapq and 10000 <= read.qlen:
            r = NanoporeRead(read.query_name, read.reference_name, read.pos, (read.pos + read.qlen), read.mapq,
                             fastq_seq=read.seq)
            r.detect_snps(snps)
            if not r.snps == []:
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


if __name__ == "__main__":
    imprinted_regions = read_imprinted_data("../data/ip_gene_pos.txt")
    all_snps = load_VCF("../data/chr13.vcf")  # 50745 SNPS
    reads = load_sam_file("../data/chr13.sam", "13", all_snps)  # 8969 filtered READS out of 50581 total rs
    # Find READS overlapping with any human imprinted region
    o = get_overlapped_reads(reads, imprinted_regions)  # 86
    save_objects("../data/chr13_snps.obj", all_snps)
    save_objects("../data/chr13_reads.obj", reads)
    save_objects("../data/chr13_reads_ir.obj", o)
