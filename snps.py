# ! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 26/03/2018
"""


###################
# Human SNPs data #
###################
class SNPs:
    """
    Human SNPs data.
    chrom, pos, id, ref, alt, qual, filter, info, gt
    (ignore all indel variants)
    How many of them locate in nanopore reads? in human imprinted regions?

    e.g.
    chr19 294525 . A C 0 PASS KM=8;KFP=0;KFF=0;MTD=bwa_freebayes,bwa_gatk,bwa_platypus,isaac_strelka GT 1|0
    """

    def __init__(self, chr, id, pos, ref, alt, gt):
        self.chrom = chr
        self.id = id
        self.pos = int(pos)
        self.ref = ref
        self.alt = alt
        self.mut = ref + alt
        self.gt = gt
        self.type = ""
        self.detect_type()
        self.reads = []

    def detect_reads(self, reads_data):
        for read in reads_data:
            if self.chrom == read.chrom and read.start <= self.pos <= read.end:
                self.reads.append(read)

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
        return "{}: {}\tREF:{}, ALT:{}\tTYPE:{}. has {} nanopore reads.".format(self.chrom, self.pos, self.ref, self.alt, self.type, len(self.reads))

    def __hash__(self):
        return hash((self.chrom, self.pos))

    def __eq__(self, other):
        """Override the default Equals behavior"""
        return self.chrom == other.chrom and self.pos == other.pos and self.mut == other.mut

    def __ne__(self, other):
        """Override the default Unequal behavior"""
        return self.chrom != other.chrom or self.pos != other.pos or self.mut != other.mut


def load_VCF(vcf_file, nanopore_reads=0):
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
                if chrom == "chr19" and a != b:  # only use heter snps
                    snp = SNPs(chrom, id, pos, ref, alt, gt)
                    #snp.detect_reads(nanopore_reads)
                    if not snp.type == "indel": #and snp.reads != []:
                        all_snps.append(snp)
        f.close()
        return all_snps
    except ValueError:
        raise RuntimeError("Not the right values to unpack.")
    except IOError:
        raise IOError("This vcf file is not available.")


def count_snps(all_snps):
    """
    :param all_snps:
    :return:
    """
    MUTATIONS = {"AG": 0, "AC": 0, "AT": 0, "CA": 0, "CT": 0, "CG": 0,
                 "TA": 0, "TC": 0, "TG": 0, "GA": 0, "GC": 0, "GT": 0}
    REF = {"A": 0, "T": 0, "G": 0, "C": 0}
    ALT = {"A": 0, "T": 0, "G": 0, "C": 0}

    for snp in all_snps:
        for mutation in MUTATIONS:
            if snp.mut == mutation:
                MUTATIONS[mutation] += 1

    for snp in all_snps:
        for r in REF:
            if snp.ref == r:
                REF[r] += 1
        for a in ALT:
            if snp.alt == a:
                ALT[a] += 1
    return MUTATIONS, REF, ALT


def get_snp_prob(all_snps, snp):
    """
    Given a SNPs, what's the frequency of happening
    its mutation (ref -> alt)? (ignore positions?)
    :param all_snps: list of SNPs objects
    :param snp: (SNPs) an instance
    :return: probability (float) calculated based on give vcf data.
    """
    counts, ref, alt = count_snps(all_snps)
    total = 0
    ref_t = 0
    alt_t = 0
    p = 0
    p_ref = 0
    p_alt = 0

    for key in counts:
        total += counts[key]
    if snp.mut in counts:
        p = counts[snp.mut] / total
    else:
        print("Wrong ref or alt value.")

    for r in ref:
        ref_t += ref[r]
    if snp.ref in ref:
        p_ref = ref[snp.ref] / ref_t
    else:
        print("Wrong ref or alt value.")

    for a in alt:
        alt_t += alt[a]
    if snp.alt in alt:
        p_alt = alt[snp.alt] / alt_t
    else:
        print("Wrong ref or alt value.")

    return p_ref * p, p_alt * p


def count_types(all_snps):
    """
    Count numbers of each type of SNPs.
    :return: lists
    """
    indels = []
    transitions = []
    transversions = []
    for snp in all_snps:
        if snp.type == "indel":
            indels.append(snp)
        elif snp.type == "transition":
            transitions.append(snp)
        else:
            transversions.append(snp)
    return indels, transitions, transversions


def get_positions(all_snps):
    pos = []
    for snp in all_snps:
        pos.append(snp.pos)
    return pos


########################
# Target READs segments #
########################
class OverlappedRead:
    """
    A READs object has 5 attributes:
    its sequencer id,
    chromosome number,
    start position on reference genome (gh38),
    end position on reference genome,
    and its SNPs objects.

    In this step, we only deal with nanopore reads that have at least one base
    overlapped with the known human imprinted gene regions.

    And to find SNPs located inside each reads using known human SNPs data.
    """

    def __init__(self, id, chr, pos1, pos2, seq):
        self.id = id
        self.chrom = chr
        self.start = int(pos1)
        self.end = int(pos2)
        self.seq = seq
        self.snps = []
        self.base = []
        self.state = ""

    def detect_snps(self, SNPs_data):
        """
        Find snps in one READs.
        :param SNPs_data: (list) class SNPs objects
        """
        for snp in SNPs_data:  # each SNPs instance, attr: chrom, id, pos, ref, alt, gt
            if snp.chrom == self.chrom and self.start <= snp.pos <= self.end:  # if position
                self.snps.append(snp)

    def get_base(self, pos):
        try:
            index = int(pos) - int(self.start)
            return self.seq[index]
        except IndexError:
            print("snp pos not right.")

    def get_bases(self):
        # TODO: make sure the bases are right. positions!
        """
        Get bases on SNP positions on READs sequence.
        :return: list of string, bases, ATCG
        """
        for snp in self.snps:
            index = int(snp.pos) - int(self.start)
            self.base.append(self.seq[index])
        return self.base

    def get_read_data(self):
        return self.id, self.snps, self.base

    def set_state(self, state):
        """Set hidden markov models for a READs."""
        self.state = state

    def __str__(self):
        return "{}:\t{}:{},{}".format(self.id, self.chrom, self.start, self.end)

    def __eq__(self, other):
        """Override the default Equals behavior"""
        return self.id == other.id and self.state == other.state

    def __ne__(self, other):
        """Override the default Unequal behavior"""
        return self.id != other.id or self.snps != other.snps or self.state != other.state


def process_all_reads(read_file, SNPs_data):
    """
    1.Nanopore reads data.
    2.SNPs data.

    :param read_file:
    :param SNPs_data:
    :return: (list) : a list of OverlappedReads instances (total 375 of them)
    """
    f = open(read_file, "r")
    all_reads = []
    for line in f:
        line = line.strip().split("\t")
        if len(line) == 8:
            id, gene, chr, info, read_pos, ref_pos, thrhld, seq = line
            ref_pos = ref_pos.strip("()").split(",")
            read = OverlappedRead(id, chr, ref_pos[0], ref_pos[1][1:], seq)
            read.detect_snps(SNPs_data)
            if read.snps == []:
                continue
            else:
                read.get_bases()
                all_reads.append(read)
    f.close()
    return all_reads


########################################################
# Manually find variants between ref and nanopore seqs #
########################################################
def read_ref_genome(fn, chr):
    """
    :param fn:
    :param chr: (int)
    :return:
    """
    chrs = ["chr1, chr19"]
    l = open(fn).readlines()
    region = []
    for index, line in enumerate(l):
        # line = line.strip().split()
        # print(line)
        if line.startswith(">chr{} ".format(chr)):
            region.append(index + 1)
        if line.startswith(">chr{} ".format(chr + 1)):  # need to change this part
            region.append(index)
    seq = ("".join(l[region[0]: region[1]])).replace("\n", "")  # 58617616 base pairs, chr19
    return seq


def compare_seqs(nano_seq, ref_seq):
    """
    :param nano_seq:
    :param ref_seq:
    :return:(dict) SNPs/same: key = position on the reads, value = (ref, alt)
    """
    snp = {}
    same = {}
    if len(nano_seq) != len(ref_seq):
        raise ValueError("Sequences don't have the same length.")
    else:
        for pos in range(len(nano_seq)):
            if nano_seq[pos] != ref_seq[pos]:
                snp[pos] = (nano_seq[pos], ref_seq[pos])
            elif nano_seq[pos] == ref_seq[pos]:
                same[pos] = (nano_seq[pos], ref_seq[pos])
        return snp, same


def find_snps_in_aln(ip_reads, ref_genome, chrom="chr19"):
    """
    :param chrom:
    :param ip_reads: (dict) imprinted nanopore reads results, id is key, ref positions is [3]
    :param ref_genome:
    :return:
    """
    nanopore_snp = {}
    read = {}
    for id in ip_reads:  # iterate through each nanopore READs sequence
        # get where the nanopore READs mapped to the reference genome
        ref_start, ref_end = (int(ip_reads[id][3][0]),
                              int(ip_reads[id][3][1][1:]))
        read_start, read_end = (int(ip_reads[id][2][0]),
                                int(ip_reads[id][2][1][1:]))

        aln_type = {}
        # get nanopore sequence and its correlated reference genome sequence to compare
        nanopore_seq = ip_reads[id][5][read_start - 1: read_end - 1]
        ref_seq = ref_genome[ref_start - 1: ref_end - 1]

        snp, same = compare_seqs(nanopore_seq, ref_seq)

        aln_type["SNPs"] = snp
        aln_type["same"] = same
        read[id] = aln_type
    nanopore_snp[chrom] = read
    return nanopore_snp