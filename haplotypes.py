import numpy as np
import math


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
            self.emission_probs[index, 0, ] = np.random.dirichlet(emission_prob)
            self.emission_probs[index, 1, ] = np.random.dirichlet(emission_prob)
            #print(self.emission_probs[index, 0, ])
        return self.emission_probs

    def read_log_likelihood(self, state, read):
        """
        Calculate the log value of probability of a snp observed in a model.
        :param state: one of the model, int 0 or 1.
        :param read: (NanoporeRead) snp
        :param snps_id: list of snp_id of the snps of the snp
        :return sum (log value of likelihood of a read base observed on a SNP position)
        """
        # if read.snps == [], will raise error (IndexError)
        read_llhd = 0
        for index, snp in enumerate(read.snps):
            try:
                t = self.get_emission()  # TODO: check if
                base_llhd = t[read.snps_id[index], state, self.observations.index(read.get_base(snp.pos))]
                read_llhd += math.log(base_llhd)  # math domain error
            except ValueError:
                continue
        return read_llhd

    def assign_reads(self, reads):
        """Calculate the parental-maternal likelihood of the reads.
        Expectation step in EM algorithm.
        
        Args:
            reads: list of reads.
        """
        pm_llhd = np.zeros((len(reads), 2), dtype=float)  # what does pm stands for?
        for idx, read in enumerate(reads):
            pm_llhd[idx, 0] = self.read_log_likelihood(0, read)
            pm_llhd[idx, 1] = self.read_log_likelihood(1, read)
            if pm_llhd[idx, 0] > pm_llhd[idx, 1]:
                self.read_assignments["m1"].append(idx)
                read.set_model(0)
            elif pm_llhd[idx, 0] < pm_llhd[idx, 1]:
                self.read_assignments["m2"].append(idx)
                read.set_model(1)
        return pm_llhd

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
        for snp_id in sr_dict.keys():
            snp = self.SNPs[snp_id]
            count = np.zeros((len(self.observations), 2))
            count[:] = pseudo_base
            for read_id in sr_dict[snp_id]:
                try:
                    read = reads[read_id]
                    if hard:
                        pm_type = np.argmax(pm_llhd[read_id, :])
                        count[self.observations.index(read.get_base(snp.pos)), pm_type] += 1
                    else:
                        count[self.observations.index(read.get_base(snp.pos)), 0] += (pm_llhd[read_id, 0] * 1)
                        count[self.observations.index(read.get_base(snp.pos)), 1] += (pm_llhd[read_id, 1] * 1)
                except ValueError:
                    continue
                for state in range(2):
                    self.emission_probs[snp_id, state, :] = (count[:, state] / np.sum(count[:, state]))

    def snps_assignments(self):
        """For a SNP, its likelihood of happened in 2 models,
        Which one is higher/more likely?"""
        for snp_id, snp in enumerate(self.emission_probs):
            alt_base = self.observations.index(self.SNPs[snp_id].alt)
            if self.emission_probs[snp_id, 0, alt_base] > self.emission_probs[snp_id, 1, alt_base]:
                self.snp_assignments["m1"].append(snp_id)
                self.SNPs[snp_id].set_model(0)
            elif self.emission_probs[snp_id, 0, alt_base] < self.emission_probs[snp_id, 1, alt_base]:
                self.snp_assignments["m2"].append(snp_id)
                self.SNPs[snp_id].set_model(1)
        return self.snp_assignments  # TODO: check if need this statement

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

    def get_snps_a(self):
        """SNPs assignments matrix."""
        return self.snp_assignments

    def get_reads_a(self):
        return self.read_assignments

    def save_model(self, fn):
        """Save markov model clustering results in VCF format."""
        with open(fn, "w") as f:
            f.write("CHR\t\tPOS\t\tREF\t\tALT\t\tGT\t\tmodel\n")
            for snp in self.SNPs:
                f.write("{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\n".format(snp.chrom, snp.pos,
                                                                    snp.ref, snp.alt,
                                                                    snp.gt, self.STATEs[snp.model]))


def snp_read_dict(reads):  # TODO: modify, pick two snps out for each read
    """Find the reads on a given snp position"""
    sr_dict = {}
    for read_id, read in enumerate(reads):
        for index, snp in enumerate(read.snps):
            snp_id = read.snps_id[index]
            if snp_id not in sr_dict.keys():
                sr_dict[snp_id] = [read_id]
            else:
                sr_dict[snp_id].append(read_id)

    #for read_id, read in enumerate(reads):
        #max_len = 0
        #for index, snp in enumerate(read.snps):
            #if len(snp.reads) > max_len:
                #max_len = len(snp.reads)
                #snp_id = read.snps_id[index]
                #snp_id2 = reads.snps_id[index-1]
                #if snp_id not in sr_dict.keys() or snp_id2 not in sr_dict.keys():
                    #sr_dict[snp_id] = [read_id]
                    #sr_dict[snp_id2] = [read_id]
                #else:
                    #sr_dict[snp_id].append(read_id)
                    #sr_dict[snp_id2].append(read_id)
    return sr_dict


def run_model(snps, reads, iter_num, hard=True):
    """Set iteration times, run a model.
    Alternative: iterate until theta converge."""
    model = HMM(snps, ["Parental", "Maternal"])
    model.init_emission_matrix()
    for _ in range(iter_num):
        sr_dict = snp_read_dict(reads)
        pm_llhd = model.assign_reads(reads)
        #print(model.emission_probs[9519,:])
        #print((np.sum(pm_llhd[:,1]),np.sum(pm_llhd[:,0])))
        model.update(reads, sr_dict, pm_llhd, hard, pseudo_base=1e-4)
    return model


def models_iterations(iter_num, snps, reads, hard=True):  # TODO: change the way to store assignment results
    """Generate a markov model multiple times."""
    m0 = run_model(snps, reads, 10, hard)
    for _ in range(iter_num):
        last_model = m0
        model = run_model(snps, reads, 10, hard)
        compare_assignments(last_model.get_reads_a(), model.get_reads_a())


def compare_assignments(o_a, a):
    """
    Compare SNPs assignments.
    :param o_a: (numpy ndarray) old assignments
    :param a: (numpy ndarray) current assignments
    :return: lists, elements changes between models
    """
    m1_to_m2 = []
    m2_to_m1 = []
    for i in o_a["m1"]:
        if i in a["m2"] and i not in a["m1"]:  # soft / hard clustering
            m1_to_m2.append(i)
    for j in o_a["m2"]:
        if j in a["m1"] and j not in a["m2"]:
            m2_to_m1.append(j)
    print("{} out of {} model1 SNPs assigned are now assigned to model2.".format(len(m1_to_m2), len(o_a["m1"])))
    print("{} out of {} model2 SNPs assigned are now assigned to model1.".format(len(m2_to_m1), len(o_a["m2"])))
    return m1_to_m2, m2_to_m1


def gold_standard(snps_data):
    """Input VCF file data."""
    pass
