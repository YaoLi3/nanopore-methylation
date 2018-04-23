import numpy as np
import math
import pandas as pd


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
        self.emission_probs = np.zeros((len(self.SNPs), len(self.STATEs), len(self.observations)))

        # Assignment results
        self.snp_assignments = np.zeros((len(self.SNPs), len(self.STATEs)))

    def init_emission_matrix(self):
        """
        Initialize emission probabilities of observing <A, T, G, C> bases on given SNP position.
        """
        for index, snp in enumerate(self.SNPs):
            emission_prob = [0.1, 0.1, 0.1, 0.1]
            emission_prob[self.observations.index(snp.ref)] = 0.4
            emission_prob[self.observations.index(snp.alt)] = 0.4
            self.emission_probs[index, 0, ] = np.random.dirichlet(emission_prob)
            self.emission_probs[index, 1, ] = np.random.dirichlet(emission_prob)
        return self.emission_probs

    def read_log_likelihood(self, state, read):
        """
        Calculate the log value of probability of a snp observed in a model.
        :param state: one of the model, int 0 or 1.
        :param read: (NanoporeRead) snp
        :param snps_id: list of snp_id of the snps of the snp
        :return log value of likelihood of the snp in given model
        """
        try:
            read_llhd = 0
            for index, snp in enumerate(read.snps):
                t = self.get_emission()  # !, ? have to
                base_llhd = t[read.snps_id[index], state, self.observations.index(read.get_base(snp.pos))]
                read_llhd += math.log(base_llhd)
            return read_llhd
        except ValueError:  # None value in the list, TODO: fix
            pass

    def assign_reads(self, reads):
        """Calculate the parental-maternal likelihood of the reads.
        Expectation step in EM algorithm.
        
        Args:
            reads: list of reads.
        """
        pm_llhd = np.zeros((len(reads), 2), dtype=float)
        assignments = np.zeros((len(reads), 2), dtype=float)
        for idx, read in enumerate(reads):
            pm_llhd[idx, 0] = self.read_log_likelihood(0, read)
            pm_llhd[idx, 1] = self.read_log_likelihood(1, read)
            if pm_llhd[idx, 0] > pm_llhd[idx, 1]:
                assignments[idx, 0] = 1
                read.set_model(0)
            elif pm_llhd[idx, 0] < pm_llhd[idx, 1]:
                assignments[idx, 1] = 1
                read.set_model(1)
        return pm_llhd, assignments

    def update(self, reads, sr_dict, pm_llhd, hard=True, pseudo_base=1e-1):
        """Update the probability matrix
        
        Args:
            reads (List): List of reads.
            sr_dict(Dictionary): Dictionary of snp snp map.
            pm_llhd (Float): Matrix of shape [Read_num,2].
            hard (Boolean): If update in hard manner.
            pseudo_base(Float): Probability normalized by adding this pseudo_base, 
                e.g. p[i] = (count[i]+pseudo_base)/(count[i]+4*pseudo_base).
        """
        try:
            for snp_id in sr_dict.keys():
                snp = self.SNPs[snp_id]
                count = np.zeros((len(self.observations), 2))
                count[:] = pseudo_base
                for read_id in sr_dict[snp_id]:
                    read = reads[read_id]
                    if hard:
                        pm_type = np.argmax(pm_llhd[read_id, :])
                        count[self.observations.index(read.get_base(snp.pos)), pm_type] += 1
                    else:
                        count[self.observations.index(read.get_base(snp.pos)), 0] += (pm_llhd[read_id, 0] * 1)
                        count[self.observations.index(read.get_base(snp.pos)), 1] += (pm_llhd[read_id, 1] * 1)
                for state in range(2):
                    self.emission_probs[snp_id, state, :] = (count[:, state] / np.sum(count[:, state]))
        except ValueError:
            pass

    def snps_assignments(self):
        """For a SNP, its likelihood of happened in 2 models,
        Which one is higher/more likely?"""
        for snp_id, snp in enumerate(self.emission_probs):
            alt_base = self.observations.index(self.SNPs[snp_id].alt)
            if self.emission_probs[snp_id, 0, alt_base] > self.emission_probs[snp_id, 1, alt_base]:
                self.snp_assignments[snp_id, 0] = 1
                self.SNPs[snp_id].set_model(0)
            elif self.emission_probs[snp_id, 0, alt_base] < self.emission_probs[snp_id, 1, alt_base]:
                self.snp_assignments[snp_id, 1] = 1
                self.SNPs[snp_id].set_model(1)

    def get_states(self):
        """Return hidden states of the model."""
        return self.STATEs

    def get_observations(self):
        """Return the observations of the model."""
        return self.observations

    def get_emission(self):  # reserve
        """Return emission probability of
        observed elements from different states."""
        return self.emission_probs

    def set_states(self, states):
        """Set states."""
        self.STATEs = states

    def set_observations(self, ob):
        """Set observations from states."""
        self.observations = ob

    def get_snps(self):
        """SNPs assignments matrix."""
        return self.snp_assignments

    def save_model(self, fn):
        """Save markov model clustering results in VCF format."""
        with open(fn, "w") as f:
            f.write("CHR\t\tPOS\t\tREF\t\tALT\t\tmodel\n")
            for snp in self.SNPs:
                f.write("{}\t\t{}\t\t{}\t\t{}\t\t{}\n".format(snp.chrom, snp.pos, snp.ref, snp.alt, self.STATEs[snp.model]))


def snp_read_dict(reads):
    """Find the reads on a given snp position"""
    sr_dict = {}
    for read_id, read in enumerate(reads):
        for index, snp in enumerate(read.snps):
            snp_id = read.snps_id[index]
            if snp_id not in sr_dict.keys():
                sr_dict[snp_id] = [read_id]
            else:
                sr_dict[snp_id].append(read_id)
    return sr_dict


def model_iterations(snps, reads, iter_num):
    """Set iteration times, run a model.
    Alternative: iterate until theta converge."""
    model = Hmm(snps, ["P", "M"])
    model.init_emission_matrix()
    total_assignments = np.zeros((len(snps), 2), dtype=float)  # total SNPs assignments, not reads
    for _ in range(iter_num):
        sr_dict = snp_read_dict(reads)
        pm_llhd, read_assign = model.assign_reads(reads)
        total_assignments += model.get_snps()
        model.snps_assignments()
        model.update(reads, sr_dict, pm_llhd, hard=True, pseudo_base=1e-4)
    return total_assignments


def compare_assignments(o_a, a):
    """
    Compare SNPs assignments.
    :param o_a: (numpy ndarray) old SNPs assignments matrix
    :param a: (numpy ndarray) current SNPs assignments matrix
    :return: lists, elements changes between models
    """
    m1_to_m2 = []
    m2_to_m1 = []
    for i in o_a[:, 0]:
        if i in a[:, 0]:
            m1_to_m2.append(i)

    for j in o_a[:, 1]:
        if j in a[:, 1]:
            m2_to_m1.append(j)

    print("{} SNPs assigned to model1 are now assigned to model2.".format(len(m1_to_m2)))
    print("{} SNPs assigned to model2 are now assigned to model1.".format(len(m2_to_m1)))
    return m1_to_m2, m2_to_m1
