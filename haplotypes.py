import numpy as np
import math
from scipy.optimize import minimize


class Hmm:
    def __init__(self, snp_data, read_data):
        self.SNPs = snp_data
        self.read = read_data
        self.ob = ["A", "T", "G", "C"]
        self.state = [0, 1]
        self.theta = np.zeros((len(self.SNPs), len(self.state), len(self.ob)))
        self.theta0 = self.theta[:0:]
        self.theta1 = self.theta[:1:]
        self.m0 = []
        self.m1 = []
        self.m0_llhd = []
        self.m1_llhd = []
        self.llhds = (self.m0, self.m1)  # likelihoods of both model in one iteration

    def init_emission_matrix(self):
        """

        :return:
        """
        for index, snp in enumerate(self.SNPs):
            emission_prob = [0.1, 0.1, 0.1, 0.1]
            emission_prob[self.ob.index(snp.ref)] = 10
            emission_prob[self.ob.index(snp.alt)] = 10
            self.theta[index:] = np.random.dirichlet(emission_prob, len(self.state))
        return self.theta

    def em(self, theta, all_reads):
        """

        :param theta:
        :param all_reads:
        :return:
        """
        def read_log_likelihood(t, state, r):
            read_llhd = 0
            for snp in r.snps:
                base_llhd = t[self.SNPs.index(snp), state, self.ob.index(r.get_base(snp.pos))]
                read_llhd += math.log2(base_llhd)
            return read_llhd

        # E-step
        for read in all_reads:
            l0 = read_log_likelihood(theta, 0, read)  # log probability of a read belongs to model 0
            l1 = read_log_likelihood(theta, 1, read)  # log probability of a read belongs to model 1
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

    def train_model(self):
        result_theta = minimize(self.em, self.init_emission_matrix(), args=self.read, method="BFGS")
        return result_theta.x
