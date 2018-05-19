import numpy as np
import math
from snps import find_most_reads_snp


class HMM:
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
        self.snp_assignments = {"m1": [], "m2": []}
        self.read_assignments = {"m1": [], "m2": []}

    def init_emission_matrix(self):
        """
        Initialize emission probabilities of observing <A, T, G, C> bases on given SNP position.
        """
        for index, snp in enumerate(self.SNPs):
            emission_prob = [10, 10, 10, 10]
            emission_prob[self.observations.index(snp.ref)] = 40
            emission_prob[self.observations.index(snp.alt)] = 40
            self.emission_probs[index, 0, ] = np.random.dirichlet(emission_prob)  # for model1
            self.emission_probs[index, 1, ] = np.random.dirichlet(emission_prob)  # for model2
        return self.emission_probs

    def read_log_likelihood(self, state, read, simulation=True):
        """
        Calculate the log value of probability of a snp observed in a model.
        :param simulation: boolean, if data is dummy data
        :param state: one of the model, int 0 or 1.
        :param read: (NanoporeRead) snp
        :param snps_id: list of snp_id of the snps of the snp
        :return sum (log value of likelihood of a read base observed on a SNP position)
        """
        read_llhd = 0

        for index, snp in enumerate(read.snps):
            try:
                t = self.get_emission()
                if simulation:
                    base_llhd = t[snp.id, state, self.observations.index(read.bases[snp.id])]
                else:
                    base_llhd = t[read.snps_id[index], state, self.observations.index(read.get_base(snp.pos))]
                read_llhd += math.log(base_llhd)  # math domain error
                return read_llhd
            except ValueError:
                continue

    def assign_reads(self, reads, simulation=True):
        """Calculate the parental-maternal likelihood of the reads.
        Expectation step in EM algorithm.
        
        Args:
            reads: list of reads.
            simulation: boolean, if data is dummy data
        """
        pm_llhd = np.zeros((len(reads), 2), dtype=float)
        self.read_assignments = {"m1": [], "m2": []}
        self.snp_assignments = {"m1": [], "m2": []}
        for idx, read in enumerate(reads):
            pm_llhd[idx, 0] = self.read_log_likelihood(0, read, simulation)
            pm_llhd[idx, 1] = self.read_log_likelihood(1, read, simulation)
            if pm_llhd[idx, 0] > pm_llhd[idx, 1]:
                self.read_assignments["m1"].append(idx)
                for snp_id in read.snps_id:
                    if snp_id not in self.snp_assignments["m1"]:
                        self.snp_assignments["m1"].append(snp_id)
            elif pm_llhd[idx, 0] < pm_llhd[idx, 1]:
                self.read_assignments["m2"].append(idx)
                for snp_id in read.snps_id:
                    if snp_id not in self.snp_assignments["m2"]:
                        self.snp_assignments["m2"].append(snp_id)
        return pm_llhd

    def update(self, reads, sr_dict, pm_llhd, hard=False, pseudo_base=1e-1, simulation=True):
        """Update the probability matrix
        
        Args:
            reads (List): List of reads.
            sr_dict(Dictionary): Dictionary of snp snp map.
            pm_llhd (Float): Matrix of shape [Read_num,2].
            hard (Boolean): If update in hard manner.
            pseudo_base(Float): Probability normalized by adding this pseudo_base, 
                e.g. p[i] = (count[i]+pseudo_base)/(count[i]+4*pseudo_base).
            simulation(Boolean): if data is dummy data
        """
        for snp_id in sr_dict.keys():
            snp = self.SNPs[snp_id]
            count = np.zeros((len(self.observations), 2))
            count[:] = pseudo_base
            for read_id in sr_dict[snp_id]:
                try:
                    read = reads[read_id]
                    if simulation:  # simulation, update counts for read base on each SNP locus
                        if hard:
                            pm_type = np.argmax(pm_llhd[read_id, :])
                            count[self.observations.index(read.bases[snp.id]), pm_type] += 1
                        else:
                            count[self.observations.index(read.bases[snp.id]), 0] += (pm_llhd[read_id, 0] * 1)  # add counts
                            count[self.observations.index(read.bases[snp.id]), 1] += (pm_llhd[read_id, 1] * 1)
                    else:  # real data
                        if hard:
                            pm_type = np.argmax(pm_llhd[read_id, :])
                            count[self.observations.index(read.get_base(snp.pos)), pm_type] += 1
                        else:
                            count[self.observations.index(read.get_base(snp.pos)), 0] += (pm_llhd[read_id, 0] * 1)
                            count[self.observations.index(read.get_base(snp.pos)), 1] += (pm_llhd[read_id, 1] * 1)
                except ValueError:
                    continue
                for state in range(len(self.STATEs)):
                    #print(count[:, state])
                    #print(np.sum(count[:, state]))
                    self.emission_probs[snp_id, state, :] = (count[:, state] / np.sum(count[:, state]))  # normalize
                    #print(self.emission_probs[snp_id, state, :])

    def get_states(self):
        """Return hidden states of the model."""
        return self.STATEs

    def get_observations(self):
        """Return the observations of the model."""
        return self.observations

    def get_emission(self):
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

    def get_reads(self):
        return self.read_assignments

    def save_model(self, fn):
        """Save markov model clustering results in VCF format."""
        pass


def snp_read_dict(reads, restrict=False):
    """Find the reads on a given snp position"""
    sr_dict = {}
    for read_id, read in enumerate(reads):
        if restrict:
            snp = find_most_reads_snp(read.snps)
            if snp is not None:
                sr_dict[snp.id] = read_id
            else:
                continue
        else:
            for index, snp in enumerate(read.snps):
                snp_id = read.snps_id[index]
                if snp_id not in sr_dict.keys():
                    sr_dict[snp_id] = [read_id]
                else:
                    sr_dict[snp_id].append(read_id)
    return sr_dict


def run_model(snps, reads, iter_num, simulation=True):
    """Set iteration times, run a model.
    Alternative: iterate until theta converge."""
    model = HMM(snps, ["Parental", "Maternal"])
    model.init_emission_matrix()
    for _ in range(iter_num):
        sr_dict = snp_read_dict(reads)
        pm_llhd = model.assign_reads(reads, simulation)
        #print(pm_llhd)
        #print(model.emission_probs[9519,:])
        #print((np.sum(pm_llhd[:,1]),np.sum(pm_llhd[:,0])))
        model.update(reads, sr_dict, pm_llhd, simulation=simulation)
    return model


def models_iterations(iter_num, snps, reads, hard=True):  # TODO: change the way to store assignment results
    """Generate a markov model multiple times."""
    m0 = run_model(snps, reads, 10, hard)
    for _ in range(iter_num):
        last_model = m0
        model = run_model(snps, reads, 10, hard)
        compare_assignments(last_model.get_reads(), model.get_reads())


def compare_assignments(o_a, a):
    """
    Compare SNPs assignments.
    :param o_a: (numpy ndarray) old assignments
    :param a: (numpy ndarray) current assignments
    :return: lists, elements changes between models
    """
    m1_to_m2 = []
    m2_to_m1 = []
    o_k = list(o_a.keys())
    a_k = list(a.keys())
    for i in o_a[o_k[0]]:
        if i in a[a_k[1]] and i not in a[a_k[0]]:  # soft / hard clustering
            m1_to_m2.append(i)
    for j in o_a[o_k[1]]:
        if j in a[a_k[0]] and j not in a[a_k[1]]:
            m2_to_m1.append(j)
    print("{} out of {} cluster 1 SNPs are now assigned to cluster 2.".format(len(m1_to_m2), len(o_a["m1"])))
    print("{} out of {} cluster 1 SNPs are now assigned to cluster 2.".format(len(m2_to_m1), len(o_a["m2"])))
    return m1_to_m2, m2_to_m1


def gold_standard(snps, snps_assignments):  # TODO: how to use gold standard from a VCF file?
    """Input VCF file data. Test data"""
    m1_A = []; m1_B = []
    for snp_id in snps_assignments["m1"]:
        if snps[snp_id].gt == "1|0":
            m1_A.append(snp_id)
        elif snps[snp_id].gt == "0|1":
            m1_B.append(snp_id)
    m2_A = []; m2_B = []
    for snp_id in snps_assignments["m2"]:
        if snps[snp_id].gt == "1|0":
            m2_A.append(snp_id)
        elif snps[snp_id].gt == "0|1":
            m2_B.append(snp_id)

    m1_num = len(snps_assignments["m1"])
    m2_num = len(snps_assignments["m2"])
    print("Model1: {}/{} SNPs, alt base on chrA; {}/{} SNPs: chrB.".format(len(m1_A), m1_num, len(m1_B), m1_num))
    print("Model2: {}/{} SNPs, alt base on chrA; {}/{} SNPs: chrB.".format(len(m2_A), m2_num, len(m2_B), m2_num))
