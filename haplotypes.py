import numpy as np
import math
import copy
from scipy.optimize import minimize


class Hmm:
    """
    Find the haplotype for each Nanopore sequencing read which has SNPs.
    """

    def __init__(self, snp_data, read_data, ob, states):
        self.SNPs = snp_data
        self.READs = read_data
        self.observations = ob
        self.STATEs = states
        self.theta = np.zeros((len(self.SNPs), len(self.STATEs), len(self.observations)))

        # Save Model data
        self.m0 = 0
        self.m1 = 0

        self.m0_reads = []  # reads of model 0
        self.m1_reads = []  # data of model 1
        self.old_m0_reads = []
        self.old_m1_reads = []
        self.m0_llhd = []  # save reads likelihood log values for one model
        self.m1_llhd = []

    def init_emission_matrix(self):
        """
        Initialize emission probabilities of observing <A, T, G, C> bases on given SNP position.
        """
        for index, snp in enumerate(self.SNPs):
            emission_prob = [0.1, 0.1, 0.1, 0.1]
            emission_prob[self.observations.index(snp.ref)] = 10
            emission_prob[self.observations.index(snp.alt)] = 10
            self.theta[index, 0, ] = np.random.dirichlet(emission_prob)
            self.theta[index, 1, ] = np.random.dirichlet(emission_prob)
        return self.theta

    def read_log_likelihood(self, t, state, r):
        """
        Calculate the log value of probability of a read observed in a model.
        :param t: theta, emission probs, 1-D numpy array
        :param state: one of the model, int 0 or 1.
        :param r: read
        :return log value of likelihood of the read in given model
        """
        read_llhd = 0
        for snp in r.snps:
            t = t.reshape((len(self.SNPs), len(self.STATEs), len(self.observations)))
            base_llhd = t[self.SNPs.index(snp), state, self.observations.index(r.get_base(snp.pos))]
            read_llhd += math.log2(base_llhd)
        return read_llhd

    def e(self, theta, all_reads):
        """
        EM algorithm to decide the haplotype of each read.
        :param theta:
        :param all_reads:
        :return: tuple, 2 model likelihood
        """
        # E-step
        # save assignments from last iteration
        self.old_m0_reads = copy.deepcopy(self.m0_reads)
        self.old_m1_reads = copy.deepcopy(self.m1_reads)
        # reset assignments
        self.m0_reads = []
        self.m1_reads = []
        self.m0_llhd = []
        self.m1_llhd = []

        for read in all_reads:
            l0 = self.read_log_likelihood(theta, 0, read)
            l1 = self.read_log_likelihood(theta, 1, read)
            if l0 > l1:
                self.m0_reads.append(read)  # save read objects
                self.m0_llhd.append(l0)  # save read likelihood log value
            elif l0 < l1:
                self.m1_reads.append(read)
                self.m1_llhd.append(l1)

    def m(self, theta, all_reads):
        # M-step (first half)
        m0 = -(sum(self.m0_llhd))  # model0 likelihood (neg log value)
        m1 = -(sum(self.m1_llhd))  # model1 likelihood
        print(m0, m1)
        llhd_array = np.array((m0, m1))
        #return self.m0 + self.m1 # try to pack them into one array # numpy array?
        #return llhd_array
        return m0

    def max_model_llhd(self, init_theta):
        """Find theta that can maximize the likelihood of a model given the current data.
        Note: minimize flatten numpy arrays into 1-d arrays."""
        result_theta = minimize(self.m, init_theta, args=self.READs, method="BFGS")
        return result_theta.x

    def if_converge(self):
        """Check if convergence happens."""
        if self.m0_reads == self.old_m0_reads and self.m1_reads == self.old_m1_reads:
            return True
        else:
            return False

    def model_iter(self):
        """Find the emission probability distribution on each SNP locus
        that maximize the likelihood of two haplotype models.
        For a read, which model it most likely came from?"""
        # while loop?
        init_theta = self.init_emission_matrix()
        self.theta = self.max_model_llhd(init_theta)

    # use matlab to make iteration (convergence process) plot, x = iteration time, y = 1.number of two groups? 2.theta?

############# basic methods ##############
    def get_states(self):
        """Return hidden states of the model."""
        return self.STATES

    def get_observs(self):
        """Return the observations of the model."""
        return self.observations

    def get_emission(self):
        """Return emission probability of
        observed elements from different states."""
        return self.theta

    def set_states(self, states):
        """Set states."""
        self.STATES = states

    def set_observs(self, ob):
        """Set observations from states."""
        self.observations = ob

    def get_data(self):
        """Return data."""
        return self.m0, self.m1

    def save(self):
        """Save markov model in a file."""
        pass

    def load(self):
        """Load model file into a Hmm object."""
        pass

    def __str__(self):
        return "This markov model. Amazing"


