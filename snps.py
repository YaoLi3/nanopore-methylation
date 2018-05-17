# ! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 26/03/2018
"""
import pysam


class SNP:
    """
    Human SNPs data.
    chr, pos, id, ref, alt, qual, filter, info, gt
    (ignore all indel variants)

    e.g.
    chr19 294525 . A C 0 PASS KM=8;KFP=0;KFF=0;MTD=bwa_freebayes,bwa_gatk,bwa_platypus,isaac_strelka GT 1|0
    """

    def __init__(self, chrom, snp_id, pos, ref, alt, gt):
        # Basic attr
        self.chr = chrom
        self.id = snp_id
        self.pos = int(pos)
        self.ref = ref
        self.alt = alt
        self.mut = ref + alt
        self.gt = gt

        # Calculated attr
        self.type = ""
        self.detect_type()
        self.close_snps = []

        # Nanopore reads attr
        self.reads = []  # reads snp_id mapped to this position

    def detect_mapped_reads(self, reads_data):
        """Find all the reads that map to this position."""
        for read in reads_data:
            if self.chr == read.chr and read.start <= self.pos <= read.end:
                self.reads.append(read)

    def detect_close_snps(self, snps):
        """Find SNPs that are close to this SNP on genome."""
        for snp in snps:
            if snp.chr == self.chr and abs(self.pos - snp.pos) < 100:  # within how many bp is close?
                self.close_snps.append(snp)

    def detect_type(self):
        """Determine the variant type."""
        TRANSITIONS = ["AG", "CT"]
        TRANSVERSIONS = ["AC", "AT", "CG", "GT"]
        if len(self.ref) > 1 or len(self.alt) > 1:
            self.type = "indel"
        elif (self.ref in TRANSITIONS[0] and self.alt in TRANSITIONS[0]) \
                or (self.ref in TRANSITIONS[1] and self.alt in TRANSITIONS[1]):
            self.type = "transition"
        else:
            for combo in TRANSVERSIONS:
                if self.ref != self.alt and self.ref in combo and self.alt in combo:
                    self.type = "transversion"

    def __str__(self):
        return "{}: {}\tREF:{}, ALT:{}\tTYPE:{}.".format(self.chr, self.pos,
                                                         self.ref, self.alt, self.type)

    def __hash__(self):
        return hash((self.chr, self.pos))

    def __eq__(self, other):
        """Override the default Equals behavior"""
        return self.chr == other.chr and self.pos == other.pos and self.mut == other.mut

    def __ne__(self, other):
        """Override the default Unequal behavior"""
        return self.chr != other.chr or self.pos != other.pos or self.mut != other.mut

    def __lt__(self, other):
        return len(self.reads) < len(other.reads)


def find_most_reads_snp(snps_list):
    """Find the SNP that has the most number of reads mapped to its position.
    Assume no "tie" situation."""
    max_reads_num = 0
    lens = []
    #best_snp = None
    #for snp in snps_list:
        #if max_reads_num < len(snp.reads):
            #max_reads_num = len(snp.reads)
            #best_snp = snp
    snps_list.sort()
    return snps_list[0]


def load_VCF(vcf_file):
    """
    Read a VCF file and return he data in it.
    :param vcf_file: (string) VCF file name
    :return: all_snps: (list) SNPs instances
    """
    try:
        f = open(vcf_file, "r")
        all_snps = []
        for line in f:
            if line.startswith("chr"):
                chrom, pos, snp_id, ref, alt, qual, fltr, \
                info, frmt, gt = line.strip().split("\t")
                a, b = gt.split("|")
                if a != b:  # only use het snps
                    snp = SNP(chrom, snp_id, pos, ref, alt, gt)
                    if not snp.type == "indel":
                        all_snps.append(snp)
        f.close()
        return all_snps
    except ValueError:
        raise RuntimeError("Not the right values to unpack.")
    except IOError:
        raise IOError("This vcf file is not available.")


def load_vcf_file(vcf_fn):  # TODO: convert file chr19.vcf into VCF format
    vf = pysam.VariantFile(vcf_fn)
    for snv in vf:
        print(snv.attr)


def map_reads(snps, reads):
    """
    Find reads mapped to each SNP position.
    :return list of snps have read mapped to
    """
    return [snp.detect_mapped_reads(reads) for snp in snps]
