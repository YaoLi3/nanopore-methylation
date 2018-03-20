#! /usr/bin/env python
"""
Find parental or maternal haplotypes of Nanopore reads.
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 08/Feb/2018
"""
import os
import sys
import glob
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
import numpy as np
import random


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
        self.emissions = []

        self.d0 = []  # data of state[0]
        self.d1 = []  # data of state[1]

        self.old_d0 = []
        self.old_d1 = []

    def symbol_to_index(self, s):
        """
        :param s: (str) DNA symbol
        :return: (int) index
        """
        try:
            return self.OBSERVATIONS.index(s)
        except ValueError:
            print("Wrong symbol.")

    def cal_snp_prob(self, snp, n=0):
        """
        Calculates a probability distribution <A, T, C, G> for each SNP position
        :return:
        """
        try:
            prob_dist = [0.1, 0.1, 0.1, 0.1]
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
        # TODO: when initialize, use all or half the reads to calculate emission probs
        """
        Initialize hmm model.
        Input data, randomly split SNPs data into two groups.
        Parameter A: transition probability.
        Parameter B: emission probability.
        Parameter pi ? from split_data function?
        """
        self.transition = np.array([[1, 0], [0, 1]])
        self.emission = self.init_emission()
        self.emissions.append(self.emission)
        self.d0, self.d1 = split_data(self.READS, 0.5)

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

    def cal_base_prob(self, read, pos):
        pass

    def cal_read_prob(self, read, state):
        # TODO: to calculate read prob, use what algorithm? or just simply multiply prob of each locus?
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
            if self.emission[list_pos, state, sym] != 0:  # there should not be zeros
                P *= self.emission[list_pos, state, sym]
        return P

    def assign_reads(self):
        # TODO: updating old assignment? but it's the same with the new assignment?
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
        self.old_d0 = self.d0[:]
        self.old_d1 = self.d1[:]

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
        # TODO: why emission = 0?
        """
        Calculate the emission matrix. Update emission prob matrix.
        save and update.
        Based on new Maternal and Parental reads, calculate the new emission probability.
        Probabilities of observing A, T, G, C on certain SNP position (loci) in each model.

        e.g. P(Ri|Dm0) = P(Dm0|Ri)P(Ri)/P(Dm0) ~ P(Dm0|Ri)P(Ri)
        P(Ri|Dm0): given the model0 data, the probability of observing base Ri on i th position
        i: pos on ref genome, here consider list_pos equal to i...
        Dm0: All reads(bases on SNP pos) in model0.
        P(Ri): prob of observing base Ri on i th pos
        P(Dm0|Ri): given the base Ri on i th pos, the probability of observing model0 data

        self.emission[snp_list_pos, 0, ] = [P(A slp|D m0), P(T|Dm0), P(G|Dm0), P(C|Dm0)]
        """

        for read in self.d0:
            for snp in read.snps:
                new_prob = [0, 0, 0, 0]
                new_prob[self.symbol_to_index(snp.ref)], \
                new_prob[self.symbol_to_index(snp.alt)] \
                    = get_snp_prob(self.get_snps(self.d0), snp)
                self.emission[self.SNPs.index(snp): 0:] = new_prob

        for read in self.d1:
            for snp in read.snps:
                new_prob = [0, 0, 0, 0]
                new_prob[self.symbol_to_index(snp.ref)], \
                new_prob[self.symbol_to_index(snp.alt)] \
                    = get_snp_prob(self.get_snps(self.d0), snp)
                self.emission[self.SNPs.index(snp): 1:] = new_prob
        self.emissions.append(self.emission)
        return self.emission

    def iteration_end(self):
        # TODO: need to fix
        """
        Decide if the iteration should stop.
        When the emission probability stops changing
        or the assignment of each read stops changing,
        stop iteration. The training of the model is complete.
        :return boolean
        """
        if self.emissions[-1].any() == self.emissions[-2].any():  # compare θ(t) and θ(t-1)
            return True
        else:
            return False

    def train_model(self):
        # TODO: finalize
        """in case of an infinite loop
        e.g.
        HmmHaplotypes.train_model(listOfReads)
        """
        self.initialize()
        while not self.iteration_end():
            self.assign_reads()
            self.cal_emission()
            print("Running")
        print("Stop iteration")

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
            P_0 = self.cal_read_prob(read, 0)
            P_1 = self.cal_read_prob(read, 1)
            if P_0 > P_1:
                return self.STATES[0]
            elif P_0 < P_1:
                return self.STATES[1]

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
