import numpy as np
import math
import copy
from scipy.optimize import minimize


class Hmm:
    """
    Find the haplotype for each Nanopore sequencing read which has SNPs.
    """

    def __init__(self, snp_data, states,ob = None):
        self.SNPs = snp_data
        self.STATEs = states
        if ob is None:
            ob = ['A','T','G','C']
        self.observations = ob
        self.theta = np.zeros((len(self.SNPs), len(self.STATEs), len(self.observations)))

        # Save Model data
        self.m0 = 0
        self.m1 = 0

        # self.m0_reads = []  # reads of model 0
        # self.m1_reads = []  # data of model 1
        # self.old_m0_reads = []
        # self.old_m1_reads = []
        # self.m0_llhd = []  # save reads likelihood log values for one model
        # self.m1_llhd = []
        # self.m0_loci = []
        # self.m1_loci = []

    def init_emission_matrix(self):
        """
        Initialize emission probabilities of observing <A, T, G, C> bases on given SNP position.
        """
        for index, snp in enumerate(self.SNPs):
            emission_prob = [0.1, 0.1, 0.1, 0.1]
            emission_prob[self.observations.index(snp.ref)] = 0.4
            emission_prob[self.observations.index(snp.alt)] = 0.4
            self.theta[index, 0,] = np.random.dirichlet(emission_prob)
            self.theta[index, 1,] = np.random.dirichlet(emission_prob)
        return self.theta

    def read_log_likelihood(self, state, r):
        """
        Calculate the log value of probability of a read observed in a model.
        :param state: one of the model, int 0 or 1.
        :param r: read
        :param snps_id: list of index of the snps of the read
        :return log value of likelihood of the read in given model
        """
        read_llhd = 0
        for index, snp in enumerate(r.snps):
            t = self.get_emission()
            base_llhd = t[r.snps_id[index], state, self.observations.index(r.get_base(snp.pos))]
#            print(base_llhd)
#            print(index)
#            print(snp)
            read_llhd += math.log(base_llhd)
        return read_llhd

    def assign_reads(self, reads):
        """Calculate the parental-maternal likelihood of the reads.
        Expectation step in EM algorithm.
        
        Args:
            reads: list of reads.
        """
        pm_llhd = np.zeros((len(reads),2),dtype = float)
        for idx,read in enumerate(reads):
            pm_llhd[idx,0] = self.read_log_likelihood(0,read)
            pm_llhd[idx,1] = self.read_log_likelihood(1,read)
        return pm_llhd

    def update(self,reads,sr_dict,pm_llhd,hard = True,pseudo_base = 1e-1):
        """Update the probability matrix
        
        Args:
            reads (List): List of reads.
            sr_dict(Dictionary): Dictionary of snp read map.
            pm_llhd (Float): Matrix of shape [Read_num,2].
            hard (Boolean): If update in hard manner.
            pseudo_base(Float): Probability normalized by adding this pseudo_base, 
                e.g. p[i] = (count[i]+pseudo_base)/(count[i]+4*pseudo_base).
        """
        for snp_id in sr_dict.keys():
            snp = self.SNPs[snp_id]
            count = np.zeros((len(self.observations),2))
            count[:]=pseudo_base
            for read_id in sr_dict[snp_id]:
                read = reads[read_id]
                if hard:
                    pm_type = np.argmax(pm_llhd[read_id,:])
                    count[self.observations.index(read.get_base(snp.pos)),pm_type]+=1
                else:
                    count[self.observations.index(read.get_base(snp.pos)),0]+=(pm_llhd[read_id,0]*1)
                    count[self.observations.index(read.get_base(snp.pos)),1]+=(pm_llhd[read_id,1]*1)
            for state in range(2):
                self.theta[snp_id,state,:] = (count[:,state] /  np.sum(count[:,state]))
                    

    def snp_read_dict(self,reads):
        """Find the reads on a given snp position"""
        sr_dict = {}
        for read_id,read in enumerate(reads):
            for index,snp in enumerate(read.snps):
                snp_id = read.snps_id[index]
                if snp_id not in sr_dict.keys():
                    sr_dict[snp_id] = [read_id]
                else:
                    sr_dict[snp_id].append(read_id)
        return sr_dict


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

    #######Function no longer needed#######
    #     def if_array_sum_one(self, theta):
    #     for (x, y, z), base_prob in np.ndenumerate(theta):
    #         if sum(theta[x, y, ]) != 1:  # index is right
    #             raise ValueError("Probability distribution does not sum to one.")

    # def if_converge(self):
    #     """Check if convergence happens."""
    #     if self.m0_reads == self.old_m0_reads and self.m1_reads == self.old_m1_reads:
    #         return True
    #     else:
    #         return False

    # def model_iter(self):
    #     """Find the emission probability distribution on each SNP locus
    #     that maximize the likelihood of two haplotype models.
    #     For a read, which model it most likely came from?"""
    #     init_theta = self.init_emission_matrix()
    #     print("init theta", init_theta[0, 0,])
    #     model0, model1, model0_loci, model1_loci = self.assign_reads(init_theta)
    #     print("init model0", len(model0), "init model1", len(model1))
    #     print(len(model0_loci), len(model1_loci))
    #     matrix_0_1 = self.maximize_likelihood_for_each_snp_pos(init_theta, 0, model0, model0_loci)
    #     new_theta = self.maximize_likelihood_for_each_snp_pos(init_theta, 1, model1, model1_loci)
    #     print(new_theta[0, 0])

    #     print("old model0 {} new model0 {} old model1 {} new model1 {}".
    #           format(len(self.m0_reads), len(self.m1_reads),
    #                  len(self.old_m0_reads), len(self.old_m1_reads)))

    #     time = 0
    #     #while not self.if_converge():
    #     while time <= 1000:
    #         current_theta = copy.deepcopy(new_theta)
    #         print("loop: current theta {}".format(current_theta[0]))
    #         model0, model1, model0_loci, model1_loci = self.assign_reads(current_theta)
    #         print(len(model0), len(model1))
    #         print(len(model0_loci), len(model1_loci))
    #         theta_0 = self.maximize_likelihood_for_each_snp_pos(current_theta, 0, model0, model0_loci)
    #         theta_1 = self.maximize_likelihood_for_each_snp_pos(current_theta, 1, model1, model1_loci)
    #         print("LOOP: old model0 {} new model0 {} old model1 {} new model1 {}".
    #               format(len(self.m0_reads), len(self.m1_reads),
    #                      len(self.old_m0_reads), len(self.old_m1_reads)))
    #         new_theta = self.theta
    #         print("loop: new theta {}".format(new_theta[0]))
    #         time += 1
    #     print("time", time)
    #     return new_theta
    ##############

