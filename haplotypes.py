import numpy as np
import math
from scipy.optimize import minimize


class Hmm:
    def __init__(self, snp_data, read_data, ob, states):
        self.SNPs = snp_data
        self.READs = read_data
        self.observations = ob
        self.STATEs = states

        self.theta = np.zeros((len(self.SNPs), len(self.STATEs), len(self.observations)))
        self.theta0 = self.theta[:0:]
        self.theta1 = self.theta[:1:]
        self.m0 = []#reads of model 0
        self.m1 = []#data of model 1
        self.old_m0 = []
        self.old_m1 = []
        self.m0_llhd = []
        self.m1_llhd = []
        self.llhds = (self.m0, self.m1)  # likelihoods of both model in one iteration

    def init_emission_matrix(self):
        """
        Initialize emission probabilities of observing <A, T, G, C> bases on given SNP position.
        """
        for index, snp in enumerate(self.SNPs):
            emission_prob = [0.1, 0.1, 0.1, 0.1]
            emission_prob[self.observations.index(snp.ref)] = 10
            emission_prob[self.observations.index(snp.alt)] = 10
            self.theta[index:] = np.random.dirichlet(emission_prob, len(self.STATEs))
        return self.theta

    def em(self, theta, all_reads):
        """
        EM algorithm to decide the haplotype of each read.
        :param theta:
        :param all_reads:
        :return:
        """
        def read_log_likelihood(t, state, r):
            """Calculate the log value of probability of a read observed in a model."""
            read_llhd = 0
            for snp in r.snps:
                base_llhd = t[self.SNPs.index(snp), state, self.observations.index(r.get_base(snp.pos))]
                read_llhd += math.log2(base_llhd)
            return read_llhd

        # E-step
        for read in all_reads:
            l0 = read_log_likelihood(theta, 0, read)
            l1 = read_log_likelihood(theta, 1, read)
            if l0 > l1:
                self.m0.append(read)
                self.m0_llhd.append(l0)
            elif l0 < l1:
                self.m1.append(read)
                self.m1_llhd.append(l1)
        # M-step (first half)
        self.m0 = -(sum(self.m0_llhd))  # model0 likelihood (neg log value)
        self.m1 = -(sum(self.m1_llhd))  # model1 likelihood
        return self.llhds

    def max_model_llhd(self, init_theta):
        """Find theta that can maximize the likelihood of a model given the current data."""
        result_theta = minimize(self.em, init_theta, args=self.READs, method="BFGS")
        return result_theta.x

    def if_convergence(self):
        pass

    def model_iter(self):
        current_theta = np.zeros((len(self.SNPs), len(self.STATEs), len(self.observations)))
        previous_theta = np.zeros((len(self.SNPs), len(self.STATEs), len(self.observations)))

        #first iteration initiation:
        init_theta = self.init_emission_matrix()
        current_theta = self.max_model_llhd(init_theta)

        while not self.m0 == self.old_m0:
            previous_theta = current_theta.clone
            current_theta = self.max_model_llhd(previous_theta)