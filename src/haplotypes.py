import numpy as np
import math
import collections  # to sort dictionaries
import logging
from src.handlefiles import save_objects, load_objects
from src.images import plot_sr_dict
import matplotlib.pyplot as plt

class HMM:
    """
    A hidden markov model.
    Find the haplotype for SNPS based on Nanopore sequencing READS.
    """

    def __init__(self, snp_data, states, ob=None):
        self.SNPs = snp_data
        self.STATEs = states
        if ob is None:
            ob = ['A', 'T', 'G', 'C']
        self.observations = ob
        self.emission_probs = np.zeros((len(self.SNPs), len(self.STATEs), len(self.observations)))

        # model results
        self.read_assignments = {"m1": [], "m2": []}
        self.alleles = {"m1": {}, "m2": {}}

    def init_emission_matrix(self):
        """
        Initialize emission probabilities of observing <A, T, G, C> bases on given SNP position.
        """
        for index, snp in enumerate(self.SNPs):
            emission_prob = [10, 10, 10, 10]
            emission_prob[self.observations.index(snp.ref)] = 40
            emission_prob[self.observations.index(snp.alt)] = 40
            self.emission_probs[index, 0, ] = np.random.dirichlet(emission_prob)  # for state1
            self.emission_probs[index, 1, ] = np.random.dirichlet(emission_prob)  # for state2
        return self.emission_probs

    def read_log_likelihood(self, state, read):
        """
        Calculate the log value of probability of a snp observed in a model.
        :param state: one of the model, int 0 or 1.
        :param read: (NanoporeRead) snp
        :return sum (log value of likelihood of a read base observed on a SNP position)

        base_llhd = t[read.snps_id[index], state, self.observations.index(read.get_base(snp.pos))]
        """
        read_llhd = 0
        for index, snp in enumerate(read.snps):
            t = self.get_emission()
            base_llhd = t[snp.id, state, self.observations.index(read.bases[snp.id])]
            if base_llhd < 0:
                logging.error("base emission prob is negative.")
            read_llhd += math.log(base_llhd)  # math domain error, base_llhd = -8.84548072603065
        return read_llhd

    def assign_reads(self, reads):
        """Calculate the parental-maternal likelihood of the READS.
        Expectation step in EM algorithm.
        SR_DICT
        Args:
            reads: list of READS.
        """
        pm_llhd = np.zeros((len(reads), 2), dtype=float)
        self.read_assignments = {"m1": [], "m2": []}
        self.alleles = {"m1": {}, "m2": {}}
        for idx, read in enumerate(reads):
            pm_llhd[idx, 0] = self.read_log_likelihood(0, read)
            pm_llhd[idx, 1] = self.read_log_likelihood(1, read)
            if pm_llhd[idx, 0] > pm_llhd[idx, 1]:
                self.read_assignments["m1"].append(idx)
                for snp_id in read.snps_id:
                    if snp_id in self.alleles["m1"]:
                        self.alleles["m1"][snp_id].append(read.bases[snp_id])
                    else:
                        self.alleles["m1"][snp_id] = [read.bases[snp_id]]
            elif pm_llhd[idx, 0] < pm_llhd[idx, 1]:
                self.read_assignments["m2"].append(idx)
                for snp_id in read.snps_id:
                    if snp_id in self.alleles["m2"]:
                        self.alleles["m2"][snp_id].append(read.bases[snp_id])
                    else:
                        self.alleles["m2"][snp_id] = [read.bases[snp_id]]
        return pm_llhd

    @staticmethod
    def random_choose(l, p):
        """must be one-d, a list or an integer works. dictionaries not."""
        s = int(len(l) * p)
        return np.random.choice(l, size=s, replace=False)

    def update(self, reads, sr_dict, pm_llhd, hard=False, pseudo_base=1e-1):
        """Update the probability matrix
        
        Args:
            reads (List): List of READS.
            sr_dict(Dictionary): Dictionary of snp snp map.
            pm_llhd (Float): Matrix of shape [Read_num,2].
            hard (Boolean): If update in hard manner.
            pseudo_base(Float): Probability normalized by adding this pseudo_base, 
                e.diffs. p[_] = (count[_]+pseudo_base)/(count[_]+4*pseudo_base).

        hard: count[self.observations.index(read.get_base(snp.pos)), pm_type] += 1
        fuzzy: count[self.observations.index(read.get_base(snp.pos)), 0] += (pm_llhd[read_id, 0] * 1)
        """
        for snp_id in sr_dict.keys():
            count = np.zeros((len(self.observations), 2))
            count[:] = pseudo_base
            for read_id in sr_dict[snp_id]:
                read = reads[read_id]
                if hard:
                    pm_type = np.argmax(pm_llhd[read_id, :])
                    count[self.observations.index(read.bases[snp_id]), pm_type] += 1  # snp_id or snp.id, same
                else:
                    count[self.observations.index(read.bases[snp_id]), 0] += np.exp(pm_llhd[read_id, 0])
                    count[self.observations.index(read.bases[snp_id]), 1] += np.exp(pm_llhd[read_id, 1])
            for state in range(len(self.STATEs)):
                self.emission_probs[snp_id, state, :] = (count[:, state] / np.sum(count[:, state]))

    def alter_update(self, percen, reads, sr_dict, pm_llhd, pseudo_base=1e-1):
        """Only consider a part of SNPS site in a read.
        soft update manner"""
        snps = self.random_choose(list(sr_dict.keys()), percen)
        for snp_id in snps:
            snp = self.SNPs[snp_id]
            count = np.zeros((len(self.observations), 2))
            count[:] = pseudo_base
            for read_id in sr_dict[snp_id]:
                read = reads[read_id]
                count[self.observations.index(read.bases[snp.id]), 0] += np.exp(pm_llhd[read_id, 0])
                count[self.observations.index(read.bases[snp.id]), 1] += np.exp(pm_llhd[read_id, 1])
            for state in range(len(self.STATEs)):
                self.emission_probs[snp_id, state, :] = (count[:, state] / np.sum(count[:, state]))

    def get_alleles(self):
        """Return a dict."""
        return self.alleles

    def get_ordered_alleles(self):
        """Return a sorted (by allele loci) collection.OrderedDict dictionary."""
        m1a = collections.OrderedDict(sorted(self.alleles["m1"].items()))
        m2a = collections.OrderedDict(sorted(self.alleles["m2"].items()))
        return {"m1": m1a, "m2": m2a}

    def get_haplotypes(self, align=False):
        """Return two strings of haplotypes generated by two states."""
        if align:
            return self.read_results()
        else:
            h1 = ""; h2 = ""
            m1 = collections.OrderedDict(sorted(self.alleles["m1"].items()))
            m2 = collections.OrderedDict(sorted(self.alleles["m2"].items()))
            for k, v in m1.items():
                h1 += v
            for k, v in m2.items():
                h2 += v
            return h1, h2

    def get_snps_pos(self):
        """Return two lists of SNPS loci."""
        p1 = sorted(list(self.alleles["m1"].keys()))
        p2 = sorted(list(self.alleles["m2"].keys()))
        diff1 = []
        diff2 = []
        if len(p1) != len(p2):
            print("The number of SNPS in two states are different.")
            for snp in p1:
                if snp not in p2:
                    diff1.append(snp)
            print("These SNP loci only appear in state 1: {}".format(diff1))
            for snp in p2:
                if snp not in p1:
                    diff2.append(snp)
            print("These SNP loci only appear in state 2: {}".format(diff2))
        elif p1 == p2:
            print("SNPS loci generated by two states are the same.")
        print(p1)
        print(p2)
        return p1, p2

    def read_results(self):
        """
        Compare haplotype strings between 2 states.
        Find difference, identify different loci and bases.
        scenario 1: the locus was covered only in one state, the other state did not emit base on this locus, reflect by -.
        scenario 2: in both states, the locus was covered with different bases.
        """
        snps_loci = list(self.alleles["m1"].keys())
        for locus in self.alleles["m2"].keys():
            if locus not in snps_loci:
                snps_loci.append(locus)
        snps_loci.sort()
        h1 = ""; h2 = ""; blank = "-"
        for pos in snps_loci:
            if pos in self.alleles["m1"]:
                h1 += self.alleles["m1"][pos]
            else:
                h1 += blank
        for pos in snps_loci:
            if pos in self.alleles["m2"]:
                h2 += self.alleles["m2"][pos]
            else:
                h2 += blank
        return h1, h2

    def get_reads(self):
        return self.read_assignments

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


def snp_read_dict(reads, limited=False, min=5):
    """Find the READS on a given snp position"""
    sr_dict = {}
    for read_id, read in enumerate(reads):
        for index, snp in enumerate(read.snps):
            snp_id = read.snps_id[index]
            if snp_id not in sr_dict.keys():
                sr_dict[snp_id] = [read_id]
            else:
                sr_dict[snp_id].append(read_id)
    if not limited:
        return sr_dict
    else:
        sr_d = {}
        for snp in sr_dict:
            if len(sr_dict[snp]) > min:
                sr_d[snp] = sr_dict[snp]
        return sr_d


def run_model(snps, reads, iter_num, hard=False, updateAll=True):
    """Set iteration times, run a model.
    Alternative: iterate until theta converge."""
    model = HMM(snps, ["Parental", "Maternal"])
    model.init_emission_matrix()
    s = np.zeros((iter_num, 2))
    for _ in range(iter_num):
        sr_dict = snp_read_dict(reads, True, min=10)
        pm_llhd = model.assign_reads(reads)
        #print(model.emission_probs[9519,:])
        #print((np.sum(pm_llhd[:,1]),np.sum(pm_llhd[:,0])))
        s[_, 0] = np.sum(pm_llhd[:, 0])
        s[_, 1] = np.sum(pm_llhd[:, 1])
        if updateAll:
            model.update(reads, sr_dict, pm_llhd, hard=hard)
        else:
            model.alter_update(reads, sr_dict, pm_llhd, 0.4)
    return model


def compare_haplotypes(h11, h12):  # sr_dict
    """
    Compare two haplotype dictionaries.
    :param h11: (dict) {snp_pos: base, snp_base: base}
    :param h12: (dict) {snp_pos: base, snp_base: base}
    :return: diff (dict) Missing, Un-match SNP sites and read bases.
            same (float) num of (h11 & h12) same / num of all h11 alleles

            save same & diff snp sites
            need: snp pos (snp id), all READS map to that snp
    """
    same = {}
    diff = {}
    for snp_id in h11:
        if snp_id in h12:
            if set(h11[snp_id]) == set(h12[snp_id]):  # if bases mapped to this pos are identical
                same[snp_id] = len(h11[snp_id])
            else:  # if bases mapped to this site are different
                diff[snp_id] = (len(h11[snp_id]), len(h12[snp_id]))  # TODO: what info should be store here
        else:  # if this SNP site only occurs in h11
            diff[snp_id] = len(h11[snp_id])
    for snp in h12:
        if snp not in h11:
            diff[snp] = len(h12[snp])
    #same = h11.items() & h12.items()
    #diff = h11.items() - h12.items() | h12.items() - h11.items()
    plot_sr_dict(same, "same.png")
    plot_sr_dict(diff, "diff.png")
    return same, diff, (len(same)/len(h11), len(same)/len(h12))


def compare_models(m1, m2):
    """
    Compare haplotype clusters generate by two hidden markov models.
    :param m1: (HMM object) a HMM hidden markov model
    :param m2: (HMM) a hidden markov model
    :return (dict) snp_pos: reads numbers
    """
    snp_sites = []
    a1 = m1.get_alleles()
    a2 = m2.get_alleles()
    for state in a1:
        for s in a2:
            print("model1: state {}; model2: state {}".format(state, s))
            r = compare_haplotypes(a1[state], a2[s])
            print("On {} SNP sites, alleles are identical. {}% of model1. {}% of model2".format(len(r[0]), r[2][0], r[2][1]))
            print("On {} SNP sites, alleles are different.".format(len(r[1])))
            snp_sites.append(r[0])

    same_snp_site = {}
    for i, snp_list in enumerate(snp_sites):
        for pos in snp_list:
            if pos in snp_sites[i-1]:
                same_snp_site[pos] = snp_list[pos]
    return same_snp_site


def models_iterations(iter_times, snps, reads, hard=False, updateAll=True):
    """Generate a markov model multiple times."""
    model = run_model(snps, reads, 10, hard=hard)
    results = np.zeros((iter_times, iter_times))   # save comparision results of multiple models
    for i in range(iter_times+1):  # i = base model
        last_model = model
        model = run_model(snps, reads, 10, hard=hard, updateAll=updateAll)
        r = compare_models(last_model, model)
        results[i, i] = 1
        results[i, i+1] = r[0]
    return results


def gold_standard(g_dict, m_dict):
    """compare model dict with this gold standard dict
    only use SNPS loci that been covered by Nanopore read in the model
    no need to use all the SNPS from the VCF file"""
    diff = {}
    for key in g_dict:
        diff[key] = {}
        for snp_pos in g_dict[key]:
            if snp_pos in m_dict["m1"]:
                diff[key]["m1"] = {}
                if m_dict["m1"][snp_pos] != g_dict[key][snp_pos]:
                    diff[key]["m1"][snp_pos] = (m_dict["m1"][snp_pos], g_dict[key][snp_pos])
            if snp_pos in m_dict["m2"]:
                diff[key]["m2"] = {}
                if m_dict["m2"][snp_pos] != g_dict[key][snp_pos]:
                    diff[key]["m2"][snp_pos] = (m_dict["m2"][snp_pos], g_dict[key][snp_pos])
    return diff


if __name__ == "__main__":
    SNPS = load_objects("../data/snps.obj")
    READS = load_objects("../data/reads.obj")
    #SR_DICT = snp_read_dict(READS)
    # gene DNMT1, chr19: 10133345 - 10231286 bp
    m1 = run_model(SNPS, READS, 100)
    m2 = run_model(SNPS, READS, 100)
    r1 = compare_models(m1, m2)
    m3 = run_model(SNPS, READS, 100)
    r2 = compare_models(m2, m3)
    r3 = compare_models(m1, m3)
    print(1 - len(set(r1)-set(r2)) / len(r1))
    print(1 - len(set(r1)-set(r3)) / len(r1))
    print(1 - len(set(r2)-set(r3)) / len(r2))
