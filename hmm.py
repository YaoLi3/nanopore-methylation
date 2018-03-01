#! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 08/02/2018
"""
import numpy as np
import random


##############
# Human SNPs #
##############
class SNPs:
    """
    Handle snps
    """

    def __init__(self):
        self.SNPs = {}

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
                    h1, h2 = gt.split("|")
                    if chrom == "chr19":
                        self.SNPs[cnt] = (chrom, int(pos), id, ref, alt,
                                      qual, filter, info, format, h1, h2)
                        cnt += 1
            f.close()
            return self.SNPs  # return dataset
        except ValueError:
            raise RuntimeError("Not the right values to unpack.")
        except IOError:
            raise IOError("This vcf file is not available.")

    def process_data(self):
        """
        Extract suitable data for hmm training.
        :param self.SNPs: (dictionary) a dataset loaded from a vcf file.
        :return: data (list) contains only position and haplotypes information.
        """
        simple_data = []
        for key in self.SNPs:
            simple_data.append((self.SNPs[key][1],  # position
                         self.SNPs[key][9],  # h1 status
                         self.SNPs[key][10]))  # h2 status
        return simple_data


#######################
# Hidden Markov Model #
#######################
class HmmHaplotypes:
    """
    Clusters Nanopore reads into two haplotype groups, based on the SNPs they have (mapped to).
    This Hidden Markov Model calculates the possibility of a reads maps to paternal or maternal
    reference chromosome with R(ref base) or N(alternative base) SNPs.
    ps. we only consider heterozygosity (1|0, 0|1, etc.)
    P(Ri | Hj) = -||(k) P(Ri,k | Hj,k)
    Iterate => update the possibility in the model.
    AN object is a whole model?
    """

    def __init__(self, k, training, testing):
        """

        :param k:
        :param training:
        :param testing:
        """
        self.order = k
        self.train = training
        self.test = testing

        self.STATES = ["Parental", "Maternal"]
        self.OBSERVATIONS = ["A", "T", "G", "C"]

        self.trans_prob = np.zeros([len(self.STATES), len(self.STATES)])
        self.emiss_prob = np.zeros([len(self.OBSERVATIONS), len(self.OBSERVATIONS)])
        self.start = np.zeros([1, 2])

        self.container = {}

    def init_hmm(self):
        """

        :return:
        """
        self.trans_prob = {"Parental": {"Parental": 0.7, "Maternal": 0.3},
                           "Maternal": {"Parental": 0.4, "Maternal": 0.6}}
        self.emiss_prob = {"Parental": {"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25},
                           "Maternal": {"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}}
        self.start = {"Parental": 0.5, "Maternal": 0.5}

    def forward_algorithm(self, obs_seq):
        """

        :param obs_seq:
        :return:
        """
        N = self.trans_prob.shape[0]
        chain_len = len(obs_seq)

        F_prob = np.zeros((N, chain_len))
        F_prob[:, 0] = self.start * self.emiss_prob[:, obs_seq[0]]

        for t in range(1, chain_len):
            for n in range(N):
                F_prob[n, t] = np.dot(F_prob[:, t - 1], (self.trans_prob[:, n])) * self.emiss_prob[n, obs_seq[t]]

        return F_prob

    def backward_algorithm(self, obs_seq):
        """

        :param obs_seq:
        :return:
        """
        N = self.trans_prob.shape[0]
        chain_len = len(obs_seq)

        B_prob = np.zeros((N, chain_len))
        B_prob[:, -1:] = 1

        for t in reversed(range(chain_len - 1)):
            for n in range(N):
                B_prob[n, t] = np.sum(B_prob[:, t + 1] * self.trans_prob[n, :] * self.emiss_prob[:, obs_seq[t + 1]])

        return B_prob

    def viterbi_algorithm(self, obs_seq):
        """

        :param obs_seq:
        :return:
        """
        N = self.trans_prob.shape[0]
        chain_len = len(obs_seq)
        prev = np.zeros((chain_len - 1, N), dtype=int)

        # DP matrix containing max likelihood of state at a given time
        V_prob = np.zeros((N, chain_len))
        V_prob[:, 0] = self.start * self.emiss_prob[:, obs_seq[0]]

        for t in range(1, chain_len):
            for n in range(N):
                seq_probs = V_prob[:, t - 1] * self.trans_prob[:, n] * self.emiss_prob[n, obs_seq[t]]
                prev[t - 1, n] = np.argmax(seq_probs)
                V_prob[n, t] = np.max(seq_probs)
        return V_prob, prev

    def train_hmm(self):
        """
        Use Baum-Welch algorithm (unsupervised). we don't know the origin of reads.
        E-step:
        M-step:
        :return:
        """

    def set_prob(self, state, prob_n):
        """
        :param state:
        :param prob_n:
        :return:
        """
        self.container[state] = prob_n

    def get_prob(self, state):
        """
        :param state:
        :return:
        """
        return self.container[state]


########
# Data #
########
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


def get_positions(file_n):
    """
    :param file_n:
    :return:
    """
    poses = {}
    f = open(file_n, "r")
    for line in f:
        line = line.strip().split("\t")
        if len(line) == 7:
            pos = line[4].strip("()").split(",")
            poses[line[0]] = (int(pos[0]), int(pos[1][1:]), line[-1])
    return poses


def find_read_snp(read_f, snp_data):
    """
    :param read_f:
    :param snp_data:
    :return:
    """
    data = {}
    info = {}
    ip_reads_regions = get_positions(read_f)
    for id in ip_reads_regions:
        l = []
        start, end, seq = ip_reads_regions[id]
        for i in range(len(snp_data)):
            snp_p = snp_data[i][0]
            if start <= snp_p <= end:
                l.append(snp_data[i])
        if l != []:
            info["SNP"] = l
            info["Seq"] = seq
            data[id] = info
    return data