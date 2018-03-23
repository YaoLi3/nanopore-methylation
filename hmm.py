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
import random
import math
from scipy.stats import bernoulli, binom
from nanoporereads import *
from IPython.display import Image
from numpy.core.umath_tests import matrix_multiply as mm


#######################
# Hidden Markov Model #
#######################
class HmmHaplotypes:
    """
    Find the haplotype for each Nanopore sequencing READs which has SNPs.
    Expectation Maximizatio (EM) Algorithm.
    """

    def __init__(self, snps, reads, states, obs, start):
        """
        :param snps: list of SNPs objects
        :param reads: list of OverlappedReads objects
        :param states: list of hidden states
        :param obs: list of emissions
        :param start: start prob of STATEs 0
        """
        self.SNPs = snps
        self.READS = reads
        self.STATES = states
        self.OBSERVATIONS = obs
        self.START = start

        self.transition = np.zeros((len(self.STATES),
                                    len(self.STATES)))
        self.emission = np.zeros((len(self.SNPs),
                                  len(self.STATES),
                                  len(self.OBSERVATIONS)))
        self.old_emission = np.zeros((len(self.SNPs),
                                  len(self.STATES),
                                  len(self.OBSERVATIONS)))

        self.d0 = []
        self.d1 = []

        self.old_d0 = []
        self.old_d1 = []

        self.old_ll = []
        self.new_ll = []
        self.LL = []

########### Tools ##############
    def symbol_to_index(self, s):
        """
        :param s: (str) DNA symbol
        :return: (int) index
        """
        try:
            return self.OBSERVATIONS.index(s)
        except ValueError:
            print("Wrong symbol.")

    def get_snps(self, reads_list):  #TODO: should move to READ class
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

############# calc probs ############
    def cal_snp_prob(self, snp):
        """
        Define a probability distribution <A, T, C, G> for a SNP locus.
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

    def base_prob(self, state, read, pos):  #same as above method?
        """
        P(Ri|Mj) = P(Mj|Ri)*P(Ri)
        P(Ri) = count(Ri)/count(all bases)

        :param state: j∈{0, 1}
        :param read: a READs ?
        :param pos:i [ pos of SNPs on ref genome]
        :return: log2 value of the probability of a base observed on pos i in model j.
        """
        P_Ri = 0

    def read_prob(self, read, state):
        """
        Calculate a READs sequence prob.
        P(Read|Mj) = Σ log2 P(Ri|Mj), (i = 1, n)

        :param read: (OverlappedRead) object, self.snps = list of SNPs objects
        :param state: (int) 0 or 1 (index of list self.STATEs)
        :return: p: (float) log2 value of probability of a (READs) occurring in one (STATEs).
        """
        p = 0
        for list_pos, snp in enumerate(read.snps):
            sym = self.symbol_to_index(read.get_bases()[list_pos])
            if self.emission[list_pos, state, sym] != 0:  # there should not be zeros
                p += math.log2(self.emission[list_pos, state, sym])
        return p

############## prepare model ################
    def init_emission(self):
        """
        Initialize emission probabilities of
        observing <A, T, G, C> bases on given SNP position.
        """
        emission = np.zeros((len(self.SNPs),
                             len(self.STATES),
                             len(self.OBSERVATIONS)))
        for index, snp in enumerate(self.SNPs):
            emission[index] = np.random.dirichlet(self.cal_snp_prob(snp),
                                                  len(self.STATES))
        return emission

    def initialize_model(self):
        """
        Initialize the markov model.
        Parameter A: transition probability.
        Parameter B: emission probability.
        Parameter pi: start probability (0.5, 0.5)
        """
        self.transition = np.array([[1, 0], [0, 1]])
        self.emission = self.init_emission()
        self.d0, self.d1 = split_data(self.READS, self.START)

    def model_llhd(self):
        """Given two sets of data, calculate model likelihood"""
        d00 = 0
        d11 = 0
        d01 = 0
        d10 = 0
        for read in self.d0:
            l0 = self.read_prob(read, 0)
            l1 = self.read_prob(read, 1)
            d00 += l0
            d01 += l1
        for read in self.d1:
            l0 = self.read_prob(read, 0)
            l1 = self.read_prob(read, 1)
            d11 += l1
            d10 += l0
        return (d00, d01), (d11, d10)

        ############ EM train model #############
    def E_assign_reads(self):
            """
            Expectation. calculate data probabilities.
            Calculate probs of each READs in all possible models.
            Assign reads to the more likely model.

            e.g.
            Read1: (SNP1, SNP3, SNP5)
            P(Read1|Mj) > P(Read1|Mj+1)
            Assign Read1 to Mj model.
            """
            self.old_d0 = self.d0[:]
            self.old_d1 = self.d1[:]
            self.d0 = []
            self.d1 = []
            self.ws_A = []
            self.ws_B = []
            self.vs_A = []
            self.vs_B = []

            for read in self.old_d0:
                LL_0 = self.read_prob(read, 0)
                LL_1 = self.read_prob(read, 1)

                denom = np.exp(LL_0) + np.exp(LL_1)
                w_A = np.exp(LL_0) / denom
                w_B = np.exp(LL_1) / denom

                ws_A.append(w_A)
                ws_B.append(w_B)
                # used for calculating emission probs
                vs_A.append(np.dot(w_A, read))  # should not be READs, should a position or ...
                vs_B.append(np.dot(w_B, read))

                # update complete log likelihood
                self.new_ll += w_A * LL_0 + w_B * LL_1

                if LL_0 > LL_1:
                    self.d0.append(read)
                    read.set_state(self.STATES[0])
                elif LL_0 < LL_1:
                    self.d1.append(read)
                    read.set_state(self.STATES[1])
            for read in self.old_d1:
                LL_0 = self.read_prob(read, 0)
                LL_1 = self.read_prob(read, 1)
                if LL_0 > LL_1:
                    self.d0.append(read)
                    read.set_state(self.STATES[0])
                elif LL_0 < LL_1:
                    self.d1.append(read)
                    read.set_state(self.STATES[1])

    def neg_loglik(self, emission, n_reads, model_llkh, states):
            """
            :param emission:
            :param n_reads: number of reads
            :param model_llkh: likelihood of a model (sum of log values of all data in the model)
            :param states: a list of possible states/distribution
            :return:
            """
            return -np.sum([binom(n_reads, emission[z]).logpmf(x) for (x, z) in zip(model_llkh, states)])

    def maximizatio_emission(self):  # TODO: THIS STEP
            """
            Maximizatio. calculate model likelihood.
            Based on given data in two groups, find a parameter (emission probs)
            that maximizes the likelihood of the data from each possible states.

            L(θj) = Σ log2 P(Read i|θ), (i = 1, n) # probably are not using this algorithm to optimize
            θj: θ of model j

            L(θj).max
            θt < θt+1? statation point...
            """
            tol = 0.01
            # M-step: update values for parameters given current distribution
            # [EQN 2]
            self.emission[:0:] = np.sum(self.vs_A, 0) / np.sum(self.vs_A)
            self.emission[:1:] = np.sum(self.vs_B, 0) / np.sum(self.vs_B)

            if np.abs(self.new_ll - self.old_ll) < tol:
                pass
            self.old_ll = self.new_ll

            likelihood0 = sum()
            m = len(self.d0)
            llhd0 = minimize(self.neg_loglik, self.emission, args=(m, xs, zs), method='tnc', options={'maxiter': 100})
            llhd1 = minimize(self.neg_loglik, self.emission, args=(m, xs, zs), method='tnc', options={'maxiter': 100})
            # thetas should be parameter probability distributions...
            # xs should be datas
            # zs?
            ll0 = sum(map(self.read_prob(), self.d0))

    ############ EM train model #############
    def E_assign_reads(self):
            """
            Expectation. calculate data probabilities.
            Calculate probs of each READs in all possible models.
            Assign reads to the more likely model.

            e.g.
            Read1: (SNP1, SNP3, SNP5)
            P(Read1|Mj) > P(Read1|Mj+1)
            Assign Read1 to Mj model.
            """
            self.old_d0 = self.d0[:]
            self.old_d1 = self.d1[:]
            self.d0 = []
            self.d1 = []
            self.ws_A = []
            self.ws_B = []
            self.vs_A = []
            self.vs_B = []

            for read in self.old_d0:
                LL_0 = self.read_prob(read, 0)
                LL_1 = self.read_prob(read, 1)

                denom = np.exp(LL_0) + np.exp(LL_1)
                w_A = np.exp(LL_0) / denom
                w_B = np.exp(LL_1) / denom

                ws_A.append(w_A)
                ws_B.append(w_B)
                # used for calculating emission probs
                vs_A.append(np.dot(w_A, read))  # should not be READs, should a position or ...
                vs_B.append(np.dot(w_B, read))

                # update complete log likelihood
                self.new_ll += w_A * LL_0 + w_B * LL_1

                if LL_0 > LL_1:
                    self.d0.append(read)
                    read.set_state(self.STATES[0])
                elif LL_0 < LL_1:
                    self.d1.append(read)
                    read.set_state(self.STATES[1])
            for read in self.old_d1:
                LL_0 = self.read_prob(read, 0)
                LL_1 = self.read_prob(read, 1)
                if LL_0 > LL_1:
                    self.d0.append(read)
                    read.set_state(self.STATES[0])
                elif LL_0 < LL_1:
                    self.d1.append(read)
                    read.set_state(self.STATES[1])

    def neg_loglik(self, emission, n_reads, model_llkh, states):
            """
            :param emission:
            :param n_reads: number of reads
            :param model_llkh: likelihood of a model (sum of log values of all data in the model)
            :param states: a list of possible states/distribution
            :return:
            """
            return -np.sum([binom(n_reads, emission[z]).logpmf(x) for (x, z) in zip(model_llkh, states)])

    def maximizatio_emission(self):  # TODO: THIS STEP
            """
            Maximizatio. calculate model likelihood.
            Based on given data in two groups, find a parameter (emission probs)
            that maximizes the likelihood of the data from each possible states.

            L(θj) = Σ log2 P(Read i|θ), (i = 1, n) # probably are not using this algorithm to optimize
            θj: θ of model j

            L(θj).max
            θt < θt+1? statation point...
            """
            tol = 0.01
            # M-step: update values for parameters given current distribution
            # [EQN 2]
            self.emission[:0:] = np.sum(self.vs_A, 0) / np.sum(self.vs_A)
            self.emission[:1:] = np.sum(self.vs_B, 0) / np.sum(self.vs_B)

            if np.abs(self.new_ll - self.old_ll) < tol:
                pass
            self.old_ll = self.new_ll

            likelihood0 = sum()
            m = len(self.d0)
            llhd0 = minimize(self.neg_loglik, self.emission, args=(m, xs, zs), method='tnc', options={'maxiter': 100})
            llhd1 = minimize(self.neg_loglik, self.emission, args=(m, xs, zs), method='tnc', options={'maxiter': 100})
            # thetas should be parameter probability distributions...
            # xs should be datas
            # zs?
            ll0 = sum(map(self.read_prob(), self.d0))

########### EM train model #############
    def em(self):
        """i = sample, j = STATEs(A/B)"""
        ll0, ll1 = self.model_llhd()
        # find emission prbs to max ll0 and ll1.
        print(ll0)
        print("----")
        print(ll1)
        return
        xs = np.array([(5, 5), (9, 1), (8, 2), (4, 6), (7, 3)]) #same as ll_m0/ll_m1 ?? no no, xs = d0/d1
        thetas = np.array([[0.6, 0.4], [0.5, 0.5]])#theta is emission prob

        tol = 0.01
        max_iter = 100

        ll_old = 0
        for i in range(max_iter):
            ws_A = []
            ws_B = []

            vs_A = []
            vs_B = []

            ll_new = 0

            # E-step: calculate probability distributions over possible completions
            for x in xs:
                # multinomial (binomial) log likelihood
                ll_A = np.sum([x * np.log(thetas[0])])
                ll_B = np.sum([x * np.log(thetas[1])])

                # [EQN 1]
                denom = np.exp(ll_A) + np.exp(ll_B)
                w_A = np.exp(ll_A) / denom
                w_B = np.exp(ll_B) / denom

                ws_A.append(w_A)
                ws_B.append(w_B)

                # used for calculating theta
                vs_A.append(np.dot(w_A, x))
                vs_B.append(np.dot(w_B, x))

                # update complete log likelihood
                ll_new += w_A * ll_A + w_B * ll_B

            # M-step: update values for parameters given current distribution
            # [EQN 2]
            thetas[0] = np.sum(vs_A, 0) / np.sum(vs_A)
            thetas[1] = np.sum(vs_B, 0) / np.sum(vs_B)
            # print distribution of z for each x and current parameter estimate

            print("Iteration: %d" % (i + 1))
            print("theta_A = %.2f, theta_B = %.2f, ll = %.2f" % (thetas[0, 0], thetas[1, 0], ll_new))

            if np.abs(ll_new - ll_old) < tol:
                break
            ll_old = ll_new

    def iteration_end(self):
        """
        When the value of parameter (emission probability) convergences,
        iteration stops.
        """
        if self.emission.any() == self.old_emission.any():  # compare θ(t) and θ(t-1)
            return True
        else:
            return False

    def train_model(self):
        self.initialize_model()
        while not self.iteration_end():
            self.em()
            print("Running")
        print("Stop iteration")

############# model usage ##################
    def viterbi(self):
            pass

    def forward(self):
            pass

    def backward(self):
            pass

    def predict(self, read, alg="viterbi"):
            """
            Find the maximum likelihood for a READs given a markov model.
            :param read: (OverlappedReads) object
            :param alg: (str) algorithm to calculate the READs prob
            :return: one of the two states, READs prob
            """
            algorithms = {"viterbi": HmmHaplotypes.viterbi,
                          "backward": HmmHaplotypes.forward,
                          "forward": HmmHaplotypes.backward}
            if alg not in algorithms:
                raise ValueError("This algorithm does not exist.")
            else:
                P_0 = self.read_prob(read, 0)
                P_1 = self.read_prob(read, 1)
                if P_0 > P_1:
                    return self.STATES[0]
                elif P_0 < P_1:
                    return self.STATES[1]

############# basic methods ##############
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

    #hmm = HmmHaplotypes(snps_data, reads_data, ["P", "M"], ["A", "T", "G", "C"], 0.5)
    #hmm.initialize_model()
    #hmm.em()
    #hmm.train_model()
    #haplotypes = hmm.get_data()