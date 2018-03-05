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
    Human SNP data.
    chrom, pos, id, ref, alt, qual, filter, info, gt
    (ignore all indel variants)
    How many of them locate in nanopore reads? in human imprinted regions?

    e.g.
    chr19 294525 . A C 0 PASS KM=8;KFP=0;KFF=0;MTD=bwa_freebayes,bwa_gatk,bwa_platypus,isaac_strelka GT 1|0
    """

    def __init__(self):
        self.SNPs = {}
        self.snp_types = {}
        self.GT = ["0|1", "1|0", "2|0", "2|1", "0|2", "1|2", "0|0", "1|1"]  # heter and homozytyes

    def load_VCF(self, file):
        """
        Read a VCF file and return he data in it.
        :param file: (string) VCF file name
        :return: self.SNPs: (dictionary) data of the variants (SNPs), key are numbers
        """
        try:
            f = open(file, "r")
            cnt = 0
            for line in f:
                if line.startswith("chr"):
                    chrom, pos, id, ref, alt, qual, \
                    filter, info, format, gt \
                        = line.strip().split("\t")
                    a, b = gt.split("|")
                    if chrom == "chr19" and a != b:  # only use heter snps
                        self.SNPs[cnt] = (chrom, int(pos), id, ref, alt,
                                          qual, filter, info, format, int(a), int(b))
                        cnt += 1
            f.close()
            return self.SNPs  # return dataset
        except ValueError:
            raise RuntimeError("Not the right values to unpack.")
        except IOError:
            raise IOError("This vcf file is not available.")

    def get_snp(self):
        """
        Extract suitable data for hmm training.
        :return: data (list) contains only position and haplotypes information.
        """
        simple_data = []
        for key in self.SNPs:
            simple_data.append((self.SNPs[key][1],  # position
                                self.SNPs[key][3],  # reference base
                                self.SNPs[key][4],  # alternative base
                                self.SNPs[key][9],  # a status
                                self.SNPs[key][10]))  # b status
        return simple_data

    def if_indel(self, ref, alt):
        """
        Decide if a snp is an indel variation.
        :param ref: (string) a reference base
        :param alt: (string) a mutated base
        :return: boolean
        """
        if len(ref) > 1 or len(alt) > 1:
            return True
        else:
            return False

    def if_transition(self, ref, alt):
        """
        Decide if a snp is a transition variation.
        :param ref: (string) a reference base
        :param alt: (string) a mutated base
        :return: boolean
        """
        TRANSITIONS = ["AG", "CT"]
        if ref == alt:
            return False
        elif ref in TRANSITIONS[0] and alt in TRANSITIONS[0]:
            return True
        elif ref in TRANSITIONS[1] and alt in TRANSITIONS[1]:
            return True
        else:
            return False

    def if_transversion(self, ref, alt):
        """
        Decide if a snp is a transversion variation.
        :param ref: (string) a reference base
        :param alt: (string) a mutated base
        :return: boolean
        """
        TRANSVERSIONS = ["AC", "AT", "CG", "GT"]
        for conbo in TRANSVERSIONS:
            if ref != alt and ref in conbo and alt in conbo:
                return True
        return False

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

    def cal_prob(self, snp):
        """
        Given a snp, what's the frequency of happening
        its mutation (ref -> alt) on its position?
        :param snp: (position on ref genome, ref, alt)
        :return: probability (float) calculated based on give vcf data.
        """
        pass


#######################
# Hidden Markov Model #
#######################
class HmmHaplotypes:
    """
    SNPs data.
    Nanopore sequencing data.
    Cluster Nanopore reads into two haplotypes (maternal and parental)
    using SNPs inside reads as tags.
    """

    def __init__(self, snps, states, obs):
        """
        :param snps:
        :param states:
        :param obs:
        """
        self.SNPs = snps
        self.STATES = states
        self.OBSERVATIONS = obs
        self.transition = {}
        self.emission = {}
        self.M = []
        self.P = []

    def init_emission(self):
        """
        Initialize emission probabilities of observing 4 bases on given SNP position.
        For each SNP, there would a reference base and an alternative base.
        observation (A, T, G, C), order specific.

        e.g.
        SNP99: ref: C. alt: T. GT: 1|0 (GT irrelevant?)
        emission[SNP99] = random.dirichlet(0, 10, 10, 0)
        """
        for snp in self.SNPs:
            #r_p = (0, 0, 0, 0) # random/artificial probability
            pos = self.SNPs[snp][2]
            ref = self.SNPs[snp][3]
            alt = self.SNPs[snp][4]
            #index1 = self.OBSERVATIONS.index(ref)
            #index2 = self.OBSERVATIONS.index(alt)
            emiss = np.random.dirichlet((len(self.emiss), (10, 0, 0, 10))).transpose()  # ? don't need

        for n in range((10, 0, 0, 10)):
            self.emission[n].loc[pos] = (emiss[0][n], emiss[1][n], emiss[3][n], emiss[4][n])

    def initialize(self):
        """
        Initialize hmm model.
        Input data, randomly split SNPs data into two groups.
        Parameter A: transition probability.
        Parameter B: emission probability.
        Parameter pi ? from split_data function?
        """
        self.M, self.P = split_data(self.SNPs, 0.5)
        self.transition = {"Parental": {"P": 1, "M": 0}, "Maternal": {"P": 0, "M": 1}}
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

    def cal_read_prob(self, read, state):
        """
        In a given model, calculate the probability of the read which has certain SNPs occurs,
        by multiplying the probability of each SNPs in this model.
        P(R) = Π (Ri|Mi)
        :param read: nanopore read
        :param state: (string) one of the model
        :return:

        e.g.
        Read1: SNP1, SNP3, SNP5
        In Maternal model:
        P(Read1|Maternal) = P(SNP1|Maternal) * P(SNP3|Maternal) * P(SNP5|Maternal)
        """
        pass

    def cal_n_assign(self, reads):
        """
        Baum-Welch algorithm (unsupervised).
        Parameters: emission probability.
        E-step: Calculate the probability of each read in two models, respectively.
        M-step: Based on results, assign the read to the model gave the bigger probability.

        e.g.
        Read1: SNP1, SNP3, SNP5
        E-step:
        In Maternal model:
        P(Read1|Maternal) = P(SNP1|Maternal) * P(SNP3|Maternal) * P(SNP5|Maternal)
        In Parental model:
        P(Read1|Parental) = P(SNP1|Parental) * P(SNP3|Parental) * P(SNP5|Parental)
        M-step:
        P(Read1|Maternal) > P(Read1|Parental)
        Assign Read1 to Maternal model.
        """
        new_M = []
        new_P = []
        for read in reads:
            P_M = self.cal_read_prob(read, "Maternal")
            P_P = self.cal_read_prob(read, "Parental")

            if P_M > P_P:
                new_M.append(read)
            elif P_M < P_P:
                new_P.append(read)
            else:
                self.random_assign(read, new_P, new_M)
        self.set_data(new_M, new_P)

    def cal_emission(self):
        """
        Based on new Maternal and Parental reads, calculate the emission probability.
        The probability of observing A, T, C, G on certain SNP position in each model.
        """
        pass

    def iteration_end(self):
        """
        When the emission probability stops changing
        or the assignment of each read stops changing,
        stop iteration. The training of the model is complete.
        """
        return False

    def train_model(self, reads):
        """in case of an infinite loop"""
        while not self.iteration_end():
            self.cal_n_assign(reads)
            self.cal_emission()

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
        self.M = m
        self.P = p


################
# Process data #
################
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


def read_find_ip_results(fn):
    """
    Read a text file where save ... data, tab delimiter ...
    :param fn:
    :return:
    """
    f = open(fn, "r")
    data = {}
    for line in f:
        line = line.strip().split("\t")
        if len(line) == 7:
            id, gene, info, read_pos, ref_pos, thrhld, seq = line
            read_pos = read_pos.strip("()").split(",")
            ref_pos = ref_pos.strip("()").split(",")
            data[id] = (gene, info, read_pos, ref_pos, thrhld, seq)
    return data


def get_positions(file_n):
    """
    :param file_n: (string) file name
    :return: (dict) key = reads id, value = (start, end) genome positions of the read
    """
    poses = {}
    f = open(file_n, "r")
    for line in f:
        line = line.strip().split("\t")
        if len(line) == 7:
            pos = line[4].strip("()").split(",")
            poses[line[0]] = (int(pos[0]), int(pos[1][1:]))
    return poses


def find_snps_in_read(read_f, snp_data):
    """
    :param read_f: (string) find_imprinted_result file name
    :param snp_data: (list) [(snp genome position, a status, b status)]
    :return: (dict) data: key = id, value = info. key = snp, seq.
    """
    data = {}
    ip_reads_regions = get_positions(read_f)
    for id in ip_reads_regions:
        snps = []
        start, end = ip_reads_regions[id]
        for snp in snp_data:
            snp_p = snp_data[snp][1]
            if start <= snp_p <= end:
                snps.append(snp_data[snp])
        if not snps == []:
            data[id] = snps
    return data


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