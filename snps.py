# ! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 26/03/2018
"""


class SNP:
    """
    Human SNPs data.
    chrom, pos, id, ref, alt, qual, filter, info, gt
    (ignore all indel variants)
    How many of them locate in nanopore reads? in human imprinted regions?

    e.g.
    chr19 294525 . A C 0 PASS KM=8;KFP=0;KFF=0;MTD=bwa_freebayes,bwa_gatk,bwa_platypus,isaac_strelka GT 1|0
    """

    def __init__(self, chr, id, pos, ref, alt, gt):
        # Basic attr
        self.chrom = chr
        self.id = id
        self.pos = int(pos)
        self.ref = ref
        self.alt = alt
        self.mut = ref + alt
        self.gt = gt

        # Calculated attr
        self.type = ""
        self.detect_type()

        # Nanopore reads attr
        self.reads = []  # reads id mapped to this position
        self.bases = []  # read bases mapped to this position

        # Markov model attr
        self.model_llkh = [0, 0]  # likelihood for a SNP occurs in each model, respectively
        self.model = None  # SNP assignment, M1 or M2

    def detect_mapped_reads(self, reads_data):
        """Find all the reads that map to this position."""
        for read in reads_data:
            if self.chrom == read.chr and read.start <= self.pos <= read.end:
                self.reads.append(read)
                self.bases.append(read.get_base(self.pos))

    def detect_type(self):
        """
        Determine a SNP type.
        """
        TRANSITIONS = ["AG", "CT"]
        TRANSVERSIONS = ["AC", "AT", "CG", "GT"]

        if len(self.ref) > 1 or len(self.alt) > 1:
            self.type = "indel"

        elif (self.ref in TRANSITIONS[0] and self.alt in TRANSITIONS[0]) \
                or (self.ref in TRANSITIONS[1] and self.alt in TRANSITIONS[1]):
            self.type = "transition"

        else:
            for conbo in TRANSVERSIONS:
                if self.ref != self.alt and self.ref in conbo and self.alt in conbo:
                    self.type = "transversion"

    def __str__(self):
        return "{}: {}\tREF:{}, ALT:{}\tTYPE:{}. has {} nanopore reads.".format(self.chrom, self.pos, self.ref,
                                                                                self.alt, self.type, len(self.reads))

    def __hash__(self):
        return hash((self.chrom, self.pos))

    def __eq__(self, other):
        """Override the default Equals behavior"""
        return self.chrom == other.chrom and self.pos == other.pos and self.mut == other.mut

    def __ne__(self, other):
        """Override the default Unequal behavior"""
        return self.chrom != other.chrom or self.pos != other.pos or self.mut != other.mut


def load_VCF(vcf_file, nanopore_reads):
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
                chrom, pos, id, ref, alt, qual, \
                filter, info, format, gt \
                    = line.strip().split("\t")
                a, b = gt.split("|")
                if a != b:  # only use het snps
                    snp = SNP(chrom, id, pos, ref, alt, gt)
                    snp.detect_mapped_reads(nanopore_reads)
                    if not snp.type == "indel":  # and snp.reads != []:
                        all_snps.append(snp)
        f.close()
        return all_snps
    except ValueError:
        raise RuntimeError("Not the right values to unpack.")
    except IOError:
        raise IOError("This vcf file is not available.")


def get_positions(all_snps):
    """
    Return SNP positions. hg38 ref genome.
    :param all_snps: (list) a list of SNP instances
    """
    pos = []
    for snp in all_snps:
        pos.append(snp.pos)
    return pos
