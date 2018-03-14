#! /usr/bin/env python
"""
Find parental or maternal haplotypes of Nanopore reads.
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 08/Feb/2018
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
        self.mut = ref+alt
        self.gt = gt
        self.type = ""
        self.detect_type()

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
        return "{}: {}\tREF:{}, ALT:{}\tTYPE:{}".format(self.chrom, self.pos, self.ref, self.alt, self.type)

    def __eq__(self, other):
        """Override the default Equals behavior"""
        return self.chrom == other.chrom and self.pos == other.pos and self.mut == other.mut

    def __ne__(self, other):
        """Override the default Unequal behavior"""
        return self.chrom != other.chrom or self.pos != other.pos or self.mut != other.mut


def load_VCF(vcf_file, count_index=False):
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


def count_snps(all_snps):
    """
    :param all_snps:
    :return:
    """
    MUTATIONS = {"AG":0, "AC":0, "AT":0, "CA":0, "CT":0, "CG":0,
                 "TA":0, "TC":0, "TG":0, "GA":0, "GC":0, "GT":0}
    REF = {"A":0, "T":0, "G":0, "C":0}
    ALT = {"A":0, "T":0, "G":0, "C":0}

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
    Given a snp, what's the frequency of happening
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

    return p_ref*p, p_alt*p

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
        self.state = ""

    def detect_snps(self, SNPs_data):
        """
        Find snps in one read.
        :param SNPs_data: (list) class SNPs objects
        """
        for snp in SNPs_data:  # each snp instance, attr: chrom, id, pos, ref, alt, gt
            if snp.chrom == self.chrom and self.start <= snp.pos <= self.end:  # if position
                self.snps.append(snp)

    def get_bases(self):
        #TODO: make sure the bases are right. positions!
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

    def set_state(self, state):
        """Set hidden markov state for a read."""
        self.state = state

    def __str__(self):
        return "{}: {},{} \tSNPs:{}".format(self.chrom, self.start, self.end, len(self.snps))

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

        self.d0 = []  # data of state[0]
        self.d1 = []  # data of state[1]

        self.old_d0 = []
        self.old_d1 = []
        self.old_emission = np.zeros((len(self.SNPs), len(self.STATES), len(self.OBSERVATIONS)))

    def symbol_to_index(self, s):
        """
        :param s: (str) DNA symbol
        :return: (int) index
        """
        try:
            return self.OBSERVATIONS.index(s)
        except ValueError:
            print("Wrong symbol.")

    def cal_snp_prob(self, snp, n = 0):
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
        except TypeError:
            raise ValueError("Wrong type of nucleotide base in SNP.")

    def init_emission(self):
        """
        Initialize emission probabilities of observing 4 bases on given SNP position.
        For each SNP, there would a reference base and an alternative base.
        observation (A, T, G, C), order specific.

        e.g.
        SNP99: ref: C. alt: T. GT: 1|0 (GT irrelevant?)
        emission[SNP99] = random.dirichlet(0, 10, 10, 0)
        """
        emission = np.zeros((len(self.SNPs), len(self.STATES), len(self.OBSERVATIONS)))
        for index, snp in enumerate(self.SNPs):
            emission[index] = np.random.dirichlet(self.cal_snp_prob(snp), len(self.STATES))
        return emission

    def initialize(self):
        #TODO: when initialize, use all or half the reads to calculate emission probs
        """
        Initialize hmm model.
        Input data, randomly split SNPs data into two groups.
        Parameter A: transition probability.
        Parameter B: emission probability.
        Parameter pi ? from split_data function?
        """
        self.d0, self.d1 = split_data(self.READS, 0.5)
        self.transition = np.array([[1, 0], [0, 1]])
        self.emission = self.init_emission()

    def random_assign(self, r, c1, c2):
        """
        Randomly assign read into one of the two models.
        :param r: read
        :param c1: model 1
        :param c2: model 2
        """
        if random.randint(0, 1) > 0.5:
            r.set_state(self.STATES[0])
            return c1.append(r)
        else:
            r.set_state(self.STATES[1])
            return c2.append(r)

    def cal_read_prob(self, read, state):
        #TODO: to calculate read prob, use what algorithm? or just simply multiply prob of each locus?
        """
        Calculate a read sequence prob.
        In a given model, calculate the probability of the read which has certain SNPs occurs,
        by multiplying the probability of each SNPs in this model.
        P(R|Model i) = Π (R d|Model i)

        :param read: (OverlappedRead) object, .snps = list of SNPs objects
        :param state: (int) 0 or 1 (index of list self.state)
        :return: P: (float) Probability of a (read) occurring in one (state).
        """
        P = 1
        for list_pos, snp in enumerate(read.snps):
            read_symbols = read.get_bases()
            sym = self.symbol_to_index(read_symbols[list_pos])
            if self.emission[list_pos, state, sym] != 0:#there should not be zeros
                P = P * self.emission[list_pos, state, sym]
        return P

    def assign(self):
        #TODO: updating old assignment? but it's the same with the new assignment?
        """
        Calculate read probability and assign it to one state.
        Baum-Welch algorithm (unsupervised).
        Parameters: emission probability.
        E-step: Calculate the probability of each read in two models, respectively.
        M-step: Based on results, assign the read to the model gave the bigger probability.

        :param: reads: a list of OverlappedReads objects

        e.g.
        Read1: (SNP1, SNP3, SNP5)
        E-step:
        In Maternal model:
        P(Read1|Maternal) = P(Read SNP1 pos base|Maternal) * P(SNP3|Maternal) * P(SNP5|Maternal)
        In Parental model:
        P(Read1|Parental) = P(SNP1|Parental) * P(SNP3|Parental) * P(SNP5|Parental)
        M-step:
        P(Read1|Maternal) > P(Read1|Parental)
        Assign Read1 to Maternal model.
        """
        self.old_d0 = self.d0
        self.old_d1 = self.d1

        self.d0 = []
        self.d1 = []
        for read in self.old_d0:
            P_0 = self.cal_read_prob(read, 0)
            P_1 = self.cal_read_prob(read, 1)
            if P_0 > P_1:
                self.d0.append(read)
                read.set_state(self.STATES[0])
            elif P_0 < P_1:
                self.d1.append(read)
                read.set_state(self.STATES[1])

        for read in self.old_d1:
            P_0 = self.cal_read_prob(read, 0)
            P_1 = self.cal_read_prob(read, 1)
            if P_0 > P_1:
                self.d0.append(read)
                read.set_state(self.STATES[0])
            elif P_0 < P_1:
                self.d1.append(read)
                read.set_state(self.STATES[1])

    def get_snps(self, reads_list):
        """
        Extract all snps in one cluster.
        :param reads_list:
        :return: snps list (unsorted?)
        """
        snps = []
        for read in reads_list:
            for snp in read.snps:
                snps.append(snp)
        return snps

    def cal_emission(self):
        #TODO: are emission probs changing(updating) as expeceted?
        """
        Calculate the emission matrix. Update emission prob matrix.
        save and update.
        Based on new Maternal and Parental reads, calculate the new emission probability.
        Probabilities of observing A, T, G, C on certain SNP position (loci) in each model.

        e.g. P(Gi|Dm0) = P(Dm0|Gi)P(Gi)/P(Dm0) ~ P(Dm0|Gi)P(Gi)
        i: pos on ref genome, here consider list_pos equal to i...
        Dm0: All reads(bases on SNP pos) in model0.
        P(Gi): prob of observing G on i th pos

        self.emission[snp_list_pos, 0, ] = [P(A slp|D m0), P(T|Dm0), P(G|Dm0), P(C|Dm0)]
        """
        #self.old_emission = self.emission

        for read in self.d0:
            for snp in read.snps:
                new_prob = [0, 0, 0, 0]
                new_prob[self.symbol_to_index(snp.ref)],\
                new_prob[self.symbol_to_index(snp.alt)] \
                    = get_snp_prob(self.get_snps(self.d0), snp)
                self.emission[self.SNPs.index(snp), 0, ] = new_prob

        for read in self.d1:
            for snp in read.snps:
                new_prob = [0, 0, 0, 0]
                new_prob[self.symbol_to_index(snp.ref)], \
                new_prob[self.symbol_to_index(snp.alt)] \
                    = get_snp_prob(self.get_snps(self.d0), snp)
                self.emission[self.SNPs.index(snp), 1, ] = new_prob

    def iteration_end(self):
        #TODO: need to fix
        """
        Decide if the iteration should stop.
        When the emission probability stops changing
        or the assignment of each read stops changing,
        stop iteration. The training of the model is complete.
        :return boolean
        """
        if self.old_emission.any() == self.emission.any():
            return True
        else:
            return False

    def train_model(self):
        #TODO: finalize
        """in case of an infinite loop
        e.g.
        HmmHaplotypes.train_model(listOfReads)
        """
        self.initialize()
        while not self.iteration_end():
            self.cal_emission()
            self.assign()
            print("Runing")
        print("Stop")

    def viterbi(self):
        pass

    def forward(self):
        pass

    def backward(self):
        pass

    def predict(self, read, alg="viterbi"):
        """
        Find the maximum likelihood for a read given a markov model.
        :param read: (OverlappedReads) object
        :param alg: (str) algorithm to calculate the read prob
        :return: one of the two states, read prob

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
        else:
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
        """Return emission probability of
        observed elements from different states."""
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
        """Save markov model in a file."""
        pass

    def load(self):
        """Load model file into a Hmm object."""
        pass

    def __str__(self):
        return "This markov model. Amazing"


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
