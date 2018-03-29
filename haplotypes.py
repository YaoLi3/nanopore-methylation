# ! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 23/03/2018
"""
import numpy as np
import math
import copy
from scipy.optimize import minimize


class Hmm:
    """
    Find the haplotype for each Nanopore sequencing read which has SNPs.
    """

    def __init__(self, snp_data, read_data, ob, states):
        # model elements
        self.SNPs = snp_data
        self.READs = read_data
        self.bases = ob
        self.models = states

        # each snp locus, each model, each base. 3d
        self.theta = np.zeros((len(self.SNPs),
                               len(self.models),
                               len(self.bases)))

        # model data/assignments
        self.m0_reads_index = []  # reads of model 0
        self.m1_reads_index = []  # data of model 1
        self.old_m0_reads_index = []
        self.old_m1_reads_index = []
        self.snps_loci = np.zeros((len(self.models), len(self.SNPs)))  # snps assigned to each model, 2d

    ########### Calculate likelihoods ######################
    def init_emission_matrix(self):
        """
        Initialize emission probabilities of observing <A, T, G, C> bases on given SNP position.
        """
        for index, snp in enumerate(self.SNPs):
            emission_prob = [1, 1, 1, 1]
            emission_prob[self.bases.index(snp.ref)] = 10
            emission_prob[self.bases.index(snp.alt)] = 10
            self.theta[index, 0,] = np.random.dirichlet(emission_prob)  # ensure model1 and 2 init value are different
            self.theta[index, 1,] = np.random.dirichlet(emission_prob)

        return self.theta

    def one_read_llhd(self, t, state, r):
        """
        Calculate the log value of probability of a read observed in a model.
        :param t: theta, emission probs, 1-D numpy array
        :param state: one of the model, int 0 or 1.
        :param r: read
        :return log value of likelihood of the read in given model
        """
        read_llhd = 0
        for snp in r.snps:
            base_llhd = t[self.SNPs.index(snp), state, self.bases.index(r.get_base(snp.pos))]
            read_llhd += math.log2(base_llhd)

        return read_llhd

    def one_locus_llhd(self, locus_theta, locus, m_reads):
        """
        locus = snp here.
        For one SNP locus, there will be multiple read bases mapped to this position.
        Given a ATGC emission probability on this locus, a likelihood of locus can be calculated.
        :param locus_theta:
        :param locus:
        :param m_reads:
        :return:
        """
        reads_on_locus = []
        locus_llhd = 0
        for read in m_reads:  # TODO: if alternative ways are better!
            if read.start <= locus.pos <= read.end:
                reads_on_locus.append(read)
        try:
            for read in reads_on_locus:  # one single nanopore read
                base_on_locus = read.get_base(locus.pos)  # find out read base on this position
                base_prob = math.log2(locus_theta[self.bases.index(base_on_locus)])  # probability of this base emit
                print("on {} locus, the likelihood of read {} has base {} is {}".format(locus.pos,
                                                                                        read.id,
                                                                                        base_on_locus,
                                                                                        base_prob))
                locus_llhd += -base_prob  # a read, contribute to the locus likelihood
        except ValueError:
            pass

        return locus_llhd

    ######################### EM algorithm ###########################
    def assign_reads(self, emission_prob):
        """
        Expectation.
        Assign reads into a model gives it the biggest likelihood.
        Given current emission probability matrix.

        :param emission_prob: theta: an emission probability matrix
        :return: new assignments: two clusters
        """
        # save previous iteration assignments
        self.old_m0_reads_index = copy.deepcopy(self.m0_reads_index)
        self.old_m1_reads_index = copy.deepcopy(self.m1_reads_index)
        # reset assignments (reads & snps)
        self.m0_reads_index = []
        self.m1_reads_index = []
        self.snps_loci = np.zeros((len(self.models), len(self.SNPs)))
        # assign reads
        for read in self.READs:
            # save read index in all reads data, not read objects
            l0 = self.one_read_llhd(emission_prob, 0, read)
            l1 = self.one_read_llhd(emission_prob, 1, read)
            if l0 > l1:
                self.m0_reads_index.append(self.READs.index(read))
                # save snps for model0 as well
                snps_index = list(map(self.SNPs.index, read.snps))
                self.snps_loci[0, snps_index] = 1
            elif l0 < l1:
                self.m1_reads_index.append(self.READs.index(read))
                snps_index = list(map(self.SNPs.index, read.snps))
                self.snps_loci[1, snps_index] = 1

        return self.m0_reads_index, self.m1_reads_index, self.snps_loci[0], self.snps_loci[1]

    def maximize_likelihood_for_each_snp_pos(self, theta, state, m_reads, m_loci):
        """
        Maximization.
        In one model, find the emission probabilities that can give the maximum likelihood of current model.
        Use scipy.optimize.minimize to find the best the emission probability on one SNP locus.
        Update the whole emission probabilities matrix.

        :param theta:
        :param state:
        :param m_reads:
        :param m_loci:
        :return: A updated emission probability matrix given current reads assignments.
        """
        global new_emission_matrix

        def con(t):
            """ensure [A,T,G,C] probs sum to one."""
            return sum(t) - 1.0

        bnds = ((0, 1), (0, 1), (0, 1), (0, 1))  # prob value must between 0 to 1

        for locus in m_loci:  # iterate through SNPs assigned in this model
            theta_locus = theta[self.SNPs.index(locus), state,]  # a atcg distribution on this position
            # optimize, find the best prob distrib on a SNP locus
            best_locus_atgc_distrib = minimize(self.one_locus_llhd,
                                               theta_locus, args=(locus, m_reads),
                                               method="SLSQP", bounds=bnds,
                                               constraints=({'type': 'eq', 'fun': con}))
            # update the whole emission prob matrix using the best local emission prob
            new_emission_matrix = self.update_emission_matrix(locus, state, best_locus_atgc_distrib.x, theta)

            print(best_locus_atgc_distrib)
            print("best atgc probs on locus {} is {}, used to be {}".format(locus.pos,
                                                                            best_locus_atgc_distrib.x,
                                                                            theta_locus))

        return new_emission_matrix

    def update_emission_matrix(self, locus, state, best_emission_on_one_locus, previous_emission_matrix):
        """
        Update the old emission probability matrix on a single SNP position. theta[x, y, ]

        :param locus:
        :param state:
        :param best_emission_on_one_locus:
        :param previous_emission_matrix:
        :return:
        """
        locus_in_matrix = self.SNPs.index(locus)
        emission_matrix = copy.deepcopy(previous_emission_matrix)
        emission_matrix[locus_in_matrix, state] = best_emission_on_one_locus
        return emission_matrix

    ################### hmm iteration #################################
    @staticmethod
    def compare_two_clusters(old_m, m):
        """
        Compare to see if elements of two clusters are identical.
        :param old_m: list of reads
        :param m: list of reads
        :return: boolean, if exactly the same
        """
        if len(old_m) != len(m):
            return False
        else:
            return old_m == m  # TODO: how to compare two list elements?

    def if_converge(self):
        """Check if convergence happens."""
        if self.compare_two_clusters(self.old_m0_reads_index, self.m0_reads_index) \
                and self.compare_two_clusters(self.old_m1_reads_index, self.m1_reads_index):
            return True
        else:
            return False

    def one_iteration(self, current_theta, m0, m0_pos, m1, m1_pos):
        """
        Same as self.model_iter()
        :param current_theta:
        :param m0:
        :param m0_pos:
        :param m1:
        :param m1_pos:
        :return:
        """
        self.maximize_likelihood_for_each_snp_pos(current_theta, 0, m0, m0_pos)  # maximize theta
        new_theta = self.maximize_likelihood_for_each_snp_pos(current_theta, 1, m1, m1_pos)
        m0_new, m1_new, m0_pos_new, m1_pos_new = self.assign_reads(new_theta)  # assign reads
        return new_theta, m0_new, m1_new, m0_pos_new, m1_pos_new

    def model_iter(self):  # TODO: mix iteration, it is not working right now 18:50, 28/03/18
        """Find the emission probability distribution on each SNP locus
        that maximize the likelihood of two haplotype models.
        For a read, which model it most likely came from?"""
        global new_theta, new_m0, new_m1, new_m0_pos, new_m1_pos

        # Init iteration
        init_theta = self.init_emission_matrix()  # init theta value
        model0, model1, model0_loci, model1_loci = self.assign_reads(init_theta)  # init assignment
        if self.theta == init_theta:
            new_theta, new_m0, new_m1, new_m0_pos, new_m1_pos = self.one_iteration(init_theta, model0, model1,
                                                                                   model0_loci, model1_loci)
        else:
            while not self.if_converge():
                new_theta, new_m0, new_m1, new_m0_pos, new_m1_pos = self.one_iteration(new_theta, new_m0, new_m1,
                                                                                       new_m0_pos, new_m1_pos)

    ############# basic methods ##############
    def get_states(self):
        """Return hidden states of the model."""
        return self.models

    def get_observs(self):
        """Return the observations of the model."""
        return self.bases

    def get_emission(self):
        """Return emission probability of
        observed elements from different states."""
        return self.theta

    def set_states(self, states):
        """Set states."""
        self.models = states

    def set_observations(self, ob):
        """Set observations from states."""
        self.bases = ob

    def get_data(self):
        """Return data."""
        return self.m0_reads_index, self.m1_reads_index

    def save(self):
        """Save markov model in a file."""
        pass

    def load(self):
        """Load model file into a Hmm object."""
        pass

    def __str__(self):
        return "This is a markov model. You've found it."
