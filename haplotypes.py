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
        self.m0_loci = []
        self.m1_loci = []

    def init_emission_matrix(self):
        """
        Initialize emission probabilities of observing <A, T, G, C> bases on given SNP position.
        """
        for index, snp in enumerate(self.SNPs):
            emission_prob = [0.1, 0.1, 0.1, 0.1]
            emission_prob[self.observations.index(snp.ref)] = 10
            emission_prob[self.observations.index(snp.alt)] = 10
            self.theta[index, 0,] = np.random.dirichlet(emission_prob)
            self.theta[index, 1,] = np.random.dirichlet(emission_prob)
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

    def log_likelihood_read(self, t, state, r):
        """
        Calculate the log value of probability of a read observed in a model.
        :param t: theta, emission probs, 1-D numpy array
        :param state: one of the model, int 0 or 1.
        :param r: read
        :return log value of likelihood of the read in given model
        """
        read_llhd = 0
        # observations = ["A", "G", "C", "T"]
        for snp in r.snps:
            base_llhd = t[self.SNPs.index(snp), state, self.observations.index(r.get_base(snp.pos))]
            read_llhd += math.log2(base_llhd)
        return read_llhd

    def assign_reads(self, emission_prob):
        # save previous iteration assignments
        self.old_m0_reads = copy.deepcopy(self.m0_reads)
        self.old_m1_reads = copy.deepcopy(self.m1_reads)
        # reset assignments
        self.m0_reads = []
        self.m1_reads = []
        self.m0_llhd = []
        self.m1_llhd = []
        self.m0_loci = []
        self.m1_loci = []
        # assign
        for read in self.READs:
            l0 = self.log_likelihood_read(emission_prob, 0, read)
            l1 = self.log_likelihood_read(emission_prob, 1, read)
            if l0 > l1:
                self.m0_reads.append(read)  # save read objects
                self.m0_llhd.append(l0)  # save read likelihood log value
                for snp in read.snps:
                    if snp not in self.m0_loci:  # save snps for model
                        self.m0_loci.append(snp)
            elif l0 < l1:
                self.m1_reads.append(read)
                self.m1_llhd.append(l1)
                for snp in read.snps:
                    if snp not in self.m1_loci:
                        self.m1_loci.append(snp)
        return self.m0_reads, self.m1_reads, self.m0_loci, self.m1_loci

    def one_locus_llhd(self, locus_theta, locussnp, m_reads):
        reads_on_locus = []
        #observation = ["A", "G", "C", "T"]
        locus_llhd = 0
        for read in m_reads:
            if read.start <= locussnp.pos <= read.end:
                reads_on_locus.append(read)
        try:
            for read in reads_on_locus:
                base_on_locus = read.get_base(locussnp.pos)
                base_prob = math.log2(locus_theta[self.observations.index(base_on_locus)])
                #print("on {} locus, the likelihood of read {} has base {} is {}"
                    #.format(locussnp.pos, read.id, base_on_locus, base_prob))
                locus_llhd += -base_prob
        except ValueError:
            pass
        return locus_llhd

    def maximize_likelihood_for_each_snp_pos(self, theta, state, m_reads, m_loci):
        """
        For only 1 state
        Model 0
        """
        # M-step (first half)
        loci_llhd = {}
        for locus in m_loci:
            theta_locus = theta[self.SNPs.index(locus), state, ]  # a atcg distribution on this position

            locus_llhd = self.one_locus_llhd(theta_locus, locus, m_reads)
            loci_llhd[locus] = locus_llhd
            best_locus_atgc_distri = minimize(self.one_locus_llhd, theta_locus, args=(locus, m_reads), method="BFGS").x
            #print("best atgc emission prob on locus {} is {}".format(locus.pos, best_locus_atgc_distri))
            #print(best_locus_atgc_distri, theta_locus) # changes did happen on at least one base
            new_emission_matrix = self.update_emission_matrix(locus, state, best_locus_atgc_distri, theta)
        return new_emission_matrix

    def update_emission_matrix(self, locus, state, best_emission_on_one_locus, previous_emission_matirx):
        """
        For only one state.
        Model 0.
        """
        locus_in_matrix = self.SNPs.index(locus)
        emission_matrix = copy.deepcopy(previous_emission_matirx)
        emission_matrix[locus_in_matrix, state] = best_emission_on_one_locus
        return emission_matrix

    def if_array_sum_one(self, theta):
        for (x, y, z), base_prob in np.ndenumerate(theta):
            if sum(theta[x, y, ]) != 1:  # index is right
                raise ValueError("Probability distribution does not sum to one.")

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
        init_theta = self.init_emission_matrix()
        print("init theta", init_theta[0, 0,])
        model0, model1, model0_loci, model1_loci = self.assign_reads(init_theta)
        print("init model0", len(model0), "init model1", len(model1))
        print(len(model0_loci), len(model1_loci))
        matrix_0_1 = self.maximize_likelihood_for_each_snp_pos(init_theta, 0, model0, model0_loci)
        new_theta = self.maximize_likelihood_for_each_snp_pos(init_theta, 1, model1, model1_loci)
        print(new_theta[0, 0])

        print("old model0 {} new model0 {} old model1 {} new model1 {}".
              format(len(self.m0_reads), len(self.m1_reads),
                     len(self.old_m0_reads), len(self.old_m1_reads)))

        time = 0
        #while not self.if_converge():
        while time <= 1000:
            current_theta = copy.deepcopy(new_theta)
            print("loop: current theta {}".format(current_theta[0]))
            model0, model1, model0_loci, model1_loci = self.assign_reads(current_theta)
            print(len(model0), len(model1))
            print(len(model0_loci), len(model1_loci))
            theta_0 = self.maximize_likelihood_for_each_snp_pos(current_theta, 0, model0, model0_loci)
            theta_1 = self.maximize_likelihood_for_each_snp_pos(current_theta, 1, model1, model1_loci)
            print("LOOP: old model0 {} new model0 {} old model1 {} new model1 {}".
                  format(len(self.m0_reads), len(self.m1_reads),
                         len(self.old_m0_reads), len(self.old_m1_reads)))
            new_theta = self.theta
            print("loop: new theta {}".format(new_theta[0]))
            time += 1
        print("time", time)
        return new_theta

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
