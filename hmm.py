#! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 08/02/2018
"""
import numpy as np
import random


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
        self.gt = gt
        self.type = ""
        self.detect_type()

    def detect_type(self):
        """
        SNP type.
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
        return "{}: {}\tREF:{}, ALT:{}\tTYPE:{}".format(self.chrom, self.pos, self.ref, self.alt, self.type)


def load_VCF(vcf_file, count_index = False):
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
                    if not snp.type == "indel":
                        all_snps.append(snp)
        f.close()
        return all_snps
    except ValueError:
        raise RuntimeError("Not the right values to unpack.")
    except IOError:
        raise IOError("This vcf file is not available.")


def cal_prob(snp, all_snps):
    """
    Given a snp, what's the frequency of happening
    its mutation (ref -> alt) on its position?
    :param snp: (SNPs) an instance
    :param all_snps: (list)
    :return: probability (float) calculated based on give vcf data.
    """
    pass


def count_types(self):
    """
    SNPs have 3 types: indel, transition and transversion.
    This method is to determine the type of SNPs and save the results in a dictionary.
    :return: (dict) key = snp type
                    value = all the information of snp from the vcf file. e.g.
                    (chrom, int(pos), id, ref, alt, qual, filter, info, format, h1, h2))
    """
    indels = []
    transitions = []
    transversions = []
    for snp in self.SNPs:
        ref = self.SNPs[snp][3]
        alt = self.SNPs[snp][4]
        if self.if_indel(ref, alt):
            indels.append(self.SNPs[snp])
        elif self.if_transition(ref, alt):
            transitions.append(self.SNPs[snp])
        elif self.if_transversion(ref, alt):
            transversions.append(self.SNPs[snp])
    self.snp_types["indels"] = indels
    self.snp_types["transitions"] = transitions
    self.snp_types["transversions"] = transversions
    return self.snp_types


########################
# Target read segments #
########################
class OverlappedRead:
    """
    A read object has 5 attributes:
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

    def detect_snps(self, SNPs_data):
        """
        Find snps in one read.
        :param SNPs_data: (list) class SNPs objects
        """
        for snp in SNPs_data:  # each snp instance, attr: chrom, id, pos, ref, alt, gt
            if snp.chrom == self.chrom and self.start <= snp.pos <= self.end:  # if position
                self.snps.append(snp)

    def get_bases(self):
        """
        Get bases on SNP positions on read sequence.
        :return: list of string, bases, ATCG
        """
        for snp in self.snps:
            index = int(snp.pos) - int(self.start)
            self.base.append(self.seq[index])
        return self.base

    def get_read_data(self):
        return self.id, self.snps, self.base

    def __str__(self):
        return "{}: {},{} \tSNPs:{}".format(self.chrom, self.start, self.end, len(self.snps))


def process_all_reads(read_file, SNPs_data):
    """
    1.Nanopore reads data.
    2.SNPs data.

    :param read_file:
    :param SNPs_data:
    :return: (list) : OverlappedReads instances (total 375 of them)
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
            if not read.snps == []:
                read.get_bases()
                all_reads.append(read)
    f.close()
    return all_reads


#######################
# Hidden Markov Model #
#######################
class HmmHaplotypes:
    """
    Find the haplotype for a sequence.
    SNPs data.
    Nanopore sequencing data.
    Cluster Nanopore reads into two haplotypes (maternal and parental)
    using SNPs inside reads as tags.
    """

    def __init__(self, snps, reads, states, obs):
        """
        :param snps: list of SNPs objects
        :param reads: list of OverlappedReads objects
        :param states: list of hidden states
        :param obs: list of emissions
        """
        self.SNPs = snps
        self.READS = reads

        self.STATES = states
        self.OBSERVATIONS = obs

        self.start = np.zeros((1, len(self.STATES)))
        self.transition = np.zeros((len(self.STATES), len(self.STATES)))
        self.emission = np.zeros((len(self.SNPs), len(self.STATES), len(self.OBSERVATIONS)))

        self.d0 = []
        self.d1 = []

        self.old_transition = np.zeros((len(self.STATES), len(self.STATES)))
        self.old_emission = np.zeros((len(self.SNPs), len(self.STATES), len(self.OBSERVATIONS)))

    def symbol_to_index(self, s):
        """

        :param s:
        :return:
        """
        return self.OBSERVATIONS.index(s)

    def cal_snp_prob(self, snp):
        """
        Calculates a probability distribution <A, T, C, G> for each SNP position
        :return:
        """
        try:
            prob_dist = [0.000001, 0.000001, 0.000001, 0.000001]
            index1 = self.OBSERVATIONS.index(snp.ref)
            index2 = self.OBSERVATIONS.index(snp.alt)
            prob_dist[index1] = 10
            prob_dist[index2] = 10
            return prob_dist
        except ValueError:
            print("Wrong nucleotide symbol.")

    def init_emission(self):
        """
        Initialize emission probabilities of observing 4 bases on given SNP position.
        For each SNP, there would a reference base and an alternative base.
        observation (A, T, G, C), order specific.

        e.g.
        SNP99: ref: C. alt: T. GT: 1|0 (GT irrelevant?)
        emission[SNP99] = random.dirichlet(0, 10, 10, 0)
        """
        for index, snp in enumerate(self.SNPs):
            self.emission[index] = np.random.dirichlet(self.cal_snp_prob(snp), len(self.STATES))


    def initialize(self):
        """
        Initialize hmm model.
        Input data, randomly split SNPs data into two groups.
        Parameter A: transition probability.
        Parameter B: emission probability.
        Parameter pi ? from split_data function?
        """
        self.d0, self.d1 = split_data(self.READS, 0.5)
        self.transition = np.array([[1, 0], [0, 1]])
        self.init_emission()

    def random_assign(self, r, c1, c2):
        """
        Randomly assign read into one of the two models.
        :param r: read
        :param c1: model 1
        :param c2: model 2
        """
        if random.randint(0, 1) > 0.5:
            return c1.append(r)
        else:
            return c2.append(r)

    def cal_read_prob(self, read, state): # 07/03/18 night, testing here. need to extract read symbols on SNP positions.
        """
        In a given model, calculate the probability of the read which has certain SNPs occurs,
        by multiplying the probability of each SNPs in this model.
        P(R) = Π (R i|Model i)
        :param read: (OverlappedRead) object, .snps = list of SNPs objects
        :param state: (string) one of the model (int) 0 or 1
        :return:

        e.g.
        Read1: SNP1, SNP3, SNP5
        In Maternal model:
        P(Read1|Maternal) = P(SNP1|Maternal) * P(SNP3|Maternal) * P(SNP5|Maternal)
        """
        P = 0
        for list_pos, snp in enumerate(read.snps):
            i = self.SNPs.index(snp)
            read_symbols = read.get_bases()
            sym = self.symbol_to_index(read_symbols[list_pos])
            P += self.emission[i, state, sym]
        return P

    def cal_n_assign(self, reads):
        """
        Baum-Welch algorithm (unsupervised).
        Parameters: emission probability.
        E-step: Calculate the probability of each read in two models, respectively.
        M-step: Based on results, assign the read to the model gave the bigger probability.

        e.g.
        Read1: (SNP1, SNP3, SNP5)
        E-step:
        In Maternal model:
        P(Read1|Maternal) = P(SNP1|Maternal) * P(SNP3|Maternal) * P(SNP5|Maternal)
        In Parental model:
        P(Read1|Parental) = P(SNP1|Parental) * P(SNP3|Parental) * P(SNP5|Parental)
        M-step:
        P(Read1|Maternal) > P(Read1|Parental)
        Assign Read1 to Maternal model.
        """
        new_d0 = []
        new_d1 = []
        for read in reads:
            P_0 = self.cal_read_prob(read, 0)
            P_1 = self.cal_read_prob(read, 1)
            if P_0 > P_1:
                new_d0.append(read)
            elif P_0 < P_1:
                new_d1.append(read)
            else:
                self.random_assign(read, new_d0, new_d1)
        self.set_data(new_d0, new_d1)

    def cal_emission(self):
        """
        Based on new Maternal and Parental reads, calculate the emission probability.
        The probability of observing A, T, C, G on certain SNP position in each model.
        """
        self.old_emission = self.emission
        self.emission = 0

    def iteration_end(self):
        """
        When the emission probability stops changing
        or the assignment of each read stops changing,
        stop iteration. The training of the model is complete.
        """
        if self.old_emission == self.emission:  # when model parameter stops changing
            return True
        else:
            return False

    def train_model(self, reads):
        """in case of an infinite loop

        e.g.
        HmmHaplotypes.train_model(listOfReads)
        """
        while not self.iteration_end():
            self.cal_n_assign(reads)
            self.cal_emission()

    def viterbi(self):
        pass

    def forward(self):
        pass

    def backward(self):
        pass

    def predict(self, read, alg = "viterbi"):
        """
        Viterbi algorithm. Find the maximum likelihood for a read.
        :param read:
        :return: one of the two states

        e.g.
        read = (SNP1, SNP2, SNP3, SNP4)
        origin = HmmHaplotypes.predict(read)
        P(read|Maternal) = 0.7
        P(read|Parental) = 0.3
        origin is "Maternal"
        """
        algorithms = {"viterbi": HmmHaplotypes.viterbi,
                      "backward": HmmHaplotypes.forward,
                      "forward": HmmHaplotypes.backward}
        if alg not in algorithms:
            raise ValueError("This algorithm does not exist.")
        pass

    def get_states(self):
        """Return hidden states of the model."""
        return self.STATES

    def get_observs(self):
        """Return the observations of the model."""
        return self.OBSERVATIONS

    def get_trans(self):
        """Return transition probability between states."""
        return self.transition

    def get_emission(self):
        """Return emission probability of observed elements from different states."""
        return self.emission

    def set_states(self, states):
        """Set states."""
        self.STATES = states

    def set_observs(self, ob):
        """Set observations from states."""
        self.OBSERVATIONS = ob

    def set_data(self, m, p):
        """Set reads belong to two models."""
        self.d0 = m
        self.d1 = p

    def get_data(self):
        """Return data."""
        return self.d0, self.d1

    def save(self):
        pass

    def load(self):
        pass

    def __str__(self):
        return "This model."


def split_data(data, ratio):
    """
    :param data: (list) full data used to split
    :param ratio: (float) between 0 and 1. Split ratio
    :return: training_data (list) training data set for building a hidden markov model.
             test (list) testing data for a built model.
    """
    size = int(len(data) * ratio)
    training_data = []
    test = list(data)
    while len(training_data) < size:
        index = random.randrange(len(test))
        training_data.append(test.pop(index))
    return training_data, test


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
    :return:(dict) snp/same: key = position on the reads, value = (ref, alt)
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
    for id in ip_reads:  # iterate through each nanopore read sequence
        # get where the nanopore read mapped to the reference genome
        ref_start, ref_end = (int(ip_reads[id][3][0]),
                              int(ip_reads[id][3][1][1:]))
        read_start, read_end = (int(ip_reads[id][2][0]),
                                int(ip_reads[id][2][1][1:]))

        aln_type = {}
        # get nanopore sequence and its correlated reference genome sequence to compare
        nanopore_seq = ip_reads[id][5][read_start - 1: read_end - 1]
        ref_seq = ref_genome[ref_start - 1: ref_end - 1]

        snp, same = compare_seqs(nanopore_seq, ref_seq)

        aln_type["snp"] = snp
        aln_type["same"] = same
        read[id] = aln_type
    nanopore_snp[chrom] = read
    return nanopore_snp