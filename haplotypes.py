import numpy as np
import math


class Hmm:
    """
    A hidden markov model.
    Find the haplotype for SNPs based on Nanopore sequencing reads.
    """

    def __init__(self, snp_data, states, ob=None):
        self.SNPs = snp_data
        self.STATEs = states
        if ob is None:
            ob = ['A', 'T', 'G', 'C']
        self.observations = ob
        self.theta = np.zeros((len(self.SNPs), len(self.STATEs), len(self.observations)))

        # Save Model data
        self.m0 = 0  # ??
        self.m1 = 0

    def init_emission_matrix(self):
        """
        Initialize emission probabilities of observing <A, T, G, C> bases on given SNP position.
        """
        for index, snp in enumerate(self.SNPs):
            emission_prob = [0.1, 0.1, 0.1, 0.1]
            emission_prob[self.observations.index(snp.ref)] = 0.4
            emission_prob[self.observations.index(snp.alt)] = 0.4
            self.theta[index, 0, ] = np.random.dirichlet(emission_prob)
            self.theta[index, 1, ] = np.random.dirichlet(emission_prob)
        return self.theta

    def read_log_likelihood(self, state, r):
        """
        Calculate the log value of probability of a read observed in a model.
        :param state: one of the model, int 0 or 1.
        :param r: read
        :param snps_id: list of index of the snps of the read
        :return log value of likelihood of the read in given model
        """
        try:
            read_llhd = 0
            for index, snp in enumerate(r.snps):
                #print(snp)
                t = self.get_emission()  # !, ? have to
                base_llhd = t[r.snps_id[index], state, self.observations.index(r.get_base(snp.pos))]  # new way to index
                #print("base {} likelihood on this snp position is: {}".format(r.get_base(snp.pos), base_llhd))
                read_llhd += math.log(base_llhd)
            #print("read likelihood of two models: {}".format(read_llhd))
            return read_llhd
        except ValueError:
            pass

    def assign_reads(self, reads):  # new way to assign, what's changed??  just use index to save
        """Calculate the parental-maternal likelihood of the reads.
        Expectation step in EM algorithm.
        
        Args:
            reads: list of reads.
        """
        pm_llhd = np.zeros((len(reads), 2), dtype=float)
        for idx, read in enumerate(reads):
            pm_llhd[idx, 0] = self.read_log_likelihood(0, read)
            pm_llhd[idx, 1] = self.read_log_likelihood(1, read)
        return pm_llhd

    def update(self, reads, sr_dict, pm_llhd, hard=True, pseudo_base=1e-1):  # new way to optimize. MAIN
        """Update the probability matrix
        
        Args:
            reads (List): List of reads.
            sr_dict(Dictionary): Dictionary of snp read map.
            pm_llhd (Float): Matrix of shape [Read_num,2].
            hard (Boolean): If update in hard manner.
            pseudo_base(Float): Probability normalized by adding this pseudo_base, 
                e.g. p[i] = (count[i]+pseudo_base)/(count[i]+4*pseudo_base).
        """
        try:
            for snp_id in sr_dict.keys():
                snp = self.SNPs[snp_id]
                count = np.zeros((len(self.observations), 2))  # i might write similar thing before
                count[:] = pseudo_base  # what is this?
                for read_id in sr_dict[snp_id]:  # read_id = read index in list reads
                    read = reads[read_id]
                    if hard:
                        pm_type = np.argmax(pm_llhd[read_id, :])  # what is this?
                        count[self.observations.index(read.get_base(snp.pos)), pm_type] += 1
                    else:  # soft
                        count[self.observations.index(read.get_base(snp.pos)), 0] += (pm_llhd[read_id, 0] * 1)
                        count[self.observations.index(read.get_base(snp.pos)), 1] += (pm_llhd[read_id, 1] * 1)
                for state in range(2):
                    self.theta[snp_id, state, :] = (count[:, state] / np.sum(count[:, state]))
        except ValueError:
            pass

    def snp_read_dict(self, reads):  # maybe can be changed, this function is already writen elsewhere
        """Find the reads on a given snp position"""
        sr_dict = {}
        for read_id, read in enumerate(reads):  # this step is fine
            #print("snps in site this read is:".format(read.snps))
            for index, snp in enumerate(read.snps):  # why? not excuted ?
                snp_id = read.snps_id[index]
                if snp_id not in sr_dict.keys():
                    sr_dict[snp_id] = [read_id]
                else:
                    sr_dict[snp_id].append(read_id)
                    #print(read_id)
        return sr_dict

    def get_states(self):
        """Return hidden states of the model."""
        return self.STATES

    def get_observs(self):
        """Return the observations of the model."""
        return self.observations

    def get_emission(self):  # reserve
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
