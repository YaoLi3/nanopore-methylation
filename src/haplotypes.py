import numpy as np
import math
import collections
import logging
from src.handlefiles import load_objects
from src.images import plot_sr_dict, save_scatter_plot_fig
from src.nanoporereads import *


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
        self.haplotypes = {"m1": {}, "m2": {}}

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
        """l must be one-d, a list or an integer. not dictionaries."""
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

    def cluster_reads(self, reads):
        """Use current emission prob matrix to cluster any input reads."""
        pm_llhd = np.zeros((len(reads), 2), dtype=float)
        read_assignments = {"m1": [], "m2": []}
        alleles = {"m1": {}, "m2": {}}
        for idx, read in enumerate(reads):
            pm_llhd[idx, 0] = self.read_log_likelihood(0, read)
            pm_llhd[idx, 1] = self.read_log_likelihood(1, read)
            if pm_llhd[idx, 0] > pm_llhd[idx, 1]:
                read_assignments["m1"].append(idx)
                for snp_id in read.snps_id:
                    if snp_id in alleles["m1"]:
                        alleles["m1"][snp_id].append(read.bases[snp_id])
                    else:
                        alleles["m1"][snp_id] = [read.bases[snp_id]]
            elif pm_llhd[idx, 0] < pm_llhd[idx, 1]:
                read_assignments["m2"].append(idx)
                for snp_id in read.snps_id:
                    if snp_id in alleles["m2"]:
                        alleles["m2"][snp_id].append(read.bases[snp_id])
                    else:
                        alleles["m2"][snp_id] = [read.bases[snp_id]]
        h = self.get_hs(alleles)
        return h

    def get_alleles(self):
        """Return a dict."""
        return self.alleles

    def get_ordered_alleles(self):
        """Return a sorted (by allele loci) collection.OrderedDict dictionary."""
        m1a = collections.OrderedDict(sorted(self.alleles["m1"].items()))
        m2a = collections.OrderedDict(sorted(self.alleles["m2"].items()))
        return {"m1": m1a, "m2": m2a}

    @staticmethod
    def most_common(lst):
        """Find the most common elements in a list."""
        return max(set(lst), key=lst.count)

    def get_haplotypes(self):
        """When there is a tie, choose base has higher emission prob.
        Or use both in turn."""
        for key in self.alleles:
            for snp_pos in self.alleles[key]:
                base = self.most_common(self.alleles[key][snp_pos])
                self.haplotypes[key][snp_pos] = (base, len(self.alleles[key][snp_pos]))
        return self.haplotypes

    @staticmethod
    def get_hs(alleles):
        """When there is a tie, choose base has higher emission prob.
        Or use both in turn."""
        haplotypes = {"m1": {}, "m2": {}}
        for key in alleles:
            for snp_pos in alleles[key]:
                base = HMM.most_common(alleles[key][snp_pos])
                haplotypes[key][snp_pos] = (base, len(alleles[key][snp_pos]))
        return haplotypes

    def align_alleles(self):
        """
        Compare haplotype strings between 2 states.
        Find difference, identify different loci and bases.
        scenario 1: the locus was covered only in one state, the other state did not emit base on this locus, reflect by -.
        scenario 2: in both states, the locus was covered with different bases.
        """
        snps_loci = list(self.haplotypes["m1"].keys())
        for locus in self.haplotypes["m2"].keys():
            if locus not in snps_loci:
                snps_loci.append(locus)
        snps_loci.sort()
        h1 = ""; h2 = ""; blank = "-"
        for pos in snps_loci:
            if pos in self.haplotypes["m1"]:
                h1 += self.haplotypes["m1"][pos][0]
            else:
                h1 += blank
        for pos in snps_loci:
            if pos in self.haplotypes["m2"]:
                h2 += self.haplotypes["m2"][pos][0]
            else:
                h2 += blank
        return h1, h2

    def get_aligned_haplotypes(self, align=True):
        """Return two strings of haplotypes generated by two states."""
        if align:
            self.get_haplotypes()
            return self.align_alleles()
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
        return sorted(list(self.alleles["m1"].keys())), sorted(list(self.alleles["m2"].keys()))

    def get_reads(self):
        """Return reads clusters."""
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
    """
    Find the READS on a given snp position
    :param reads: list of NanoporeRead objects
    :param limited: (bool) if filter snps based on their reads num
    :param min: (int) the bare minimum number of reads mapped to a snp
    :return: (dict) {snp: read_id1, read_id2}
    """
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
            if len(sr_dict[snp]) > min:  # only SNP has more than min numbers reads mapped would be used
                sr_d[snp] = sr_dict[snp]
        return sr_d


def run_model(snps, reads, iter_num, hard=False, updateAll=True, p=0.5):
    """Set iteration times, run a model.
    Alternative: iterate until theta converge."""
    model = HMM(snps, ["Parental", "Maternal"])
    model.init_emission_matrix()
    s = np.zeros((iter_num, 2))
    for _ in range(iter_num):
        sr_dict = snp_read_dict(reads, True, min=10)
        pm_llhd = model.assign_reads(reads)#[0]
        #print(model.emission_probs[9519,:])
        #print((np.sum(pm_llhd[:,1]),np.sum(pm_llhd[:,0])))
        s[_, 0] = np.sum(pm_llhd[:, 0])
        s[_, 1] = np.sum(pm_llhd[:, 1])
        if updateAll:
            model.update(reads, sr_dict, pm_llhd, hard=hard)
        else:  # randomly choose p% SNP sites to update
            model.alter_update(reads, sr_dict, pm_llhd, p)
    #save_scatter_plot_fig(s, "likelihood_value.png")
    return model


def same_snp_sites(h1, h2):  # sr_dict
    """
    Find SNP sites with same alleles.
    """
    same = {}
    for snp_id in h1:
        if snp_id in h2 and h2[snp_id][0] == h1[snp_id][0]:
            same[snp_id] = (h1[snp_id][1], h2[snp_id][1], h1[snp_id])
    if len(h1) == 0:
        same_p = 0
    else:
        same_p = len(same) / len(h1)
    if len(h2) == 0:
        same_p2 = 0
    else:
        same_p2 = len(same) / len(h2)
    print("{}% of the first haplotype alleles are the same as the second haplotype alleles.".format(same_p))
    print("{}% of the second haplotype alleles are the same as the first haplotype alleles.".format(same_p2))
    return same


def find_common_poses(m1, m4):
    m1h = m1.get_haplotypes()
    m4h = m4.get_haplotypes()

    snp_posesm1 = []
    for snp_pos in m1h["m1"]:
        if m1h["m1"][snp_pos][1] == 1:
            pass
        else:
            if snp_pos not in snp_posesm1:
                snp_posesm1.append(snp_pos)
    snp_posesm2 = []
    for snp_pos in m1h["m2"]:
        if m1h["m2"][snp_pos][1] == 1:
            pass
        else:
            if snp_pos not in snp_posesm2:
                snp_posesm2.append(snp_pos)

    snp_poses2m1 = []
    for snp_pos in m4h["m1"]:
        if m4h["m1"][snp_pos][1] == 1:
            pass
        else:
            if snp_pos not in snp_poses2m1:
                snp_poses2m1.append(snp_pos)
    snp_poses2m2 = []
    for snp_pos in m4h["m2"]:
        if m4h["m2"][snp_pos][1] == 1:
            pass
        else:
            if snp_pos not in snp_poses2m2:
                snp_poses2m2.append(snp_pos)

    positions = [snp_posesm1, snp_posesm2, snp_poses2m1, snp_poses2m2]
    same_pos = set(positions[0])
    for s in positions[1:]:
        same_pos.intersection_update(s)

    print(len(same_pos)/len(snp_posesm1))
    print(len(same_pos)/len(snp_posesm2))
    print(len(same_pos)/len(snp_poses2m1))
    print(len(same_pos)/len(snp_poses2m2))
    print(len(same_pos))
    return same_pos


def diff_snp_sites(h1, h2):
    diff = {}
    for snp_id in h1:
        if snp_id not in h2:
            diff[snp_id] = (h1[snp_id][1], 0)
        elif h2[snp_id][0] != h1[snp_id][0]:
            diff[snp_id] = (h1[snp_id][1], h2[snp_id][1])
    for snp in h2:
        if snp not in h1:
            diff[snp] = (0, h2[snp][1])

    if len(h1) == 0:
        diff_p = 0
    else:
        diff_p = len(diff) / len(h1)
    if len(h2) == 0:
        diff_p2 = 0
    else:
        diff_p2 = len(diff) / len(h2)
    print("{}% of the first haplotype alleles are different from the second haplotype alleles.".format(diff_p))
    print("{}% of the second haplotype alleles are different from the first haplotype alleles.".format(diff_p2))
    return diff


def compare_models(a1, a2):
    """
    same. pos. read num. percentage.
    diff.
    """
    #a1 = m1.get_haplotypes()
    #a2 = m2.get_haplotypes()
    same_snp_sites(a1["m1"], a2["m1"])
    same_snp_sites(a1["m1"], a2["m2"])
    same_snp_sites(a1["m2"], a2["m1"])
    same_snp_sites(a1["m2"], a2["m1"])
    diff_snp_sites(a1["m1"], a2["m1"])
    diff_snp_sites(a1["m1"], a2["m2"])
    diff_snp_sites(a1["m2"], a2["m1"])
    diff_snp_sites(a1["m2"], a2["m2"])


def compare_clustering(m1, m2, data):
    mm1 = m1.cluster_reads(data)
    mm4 = m2.cluster_reads(data)
    same_snp_sites(mm1["m1"], mm4["m1"])
    same_snp_sites(mm1["m1"], mm4["m2"])
    same_snp_sites(mm1["m2"], mm4["m1"])
    same_snp_sites(mm1["m2"], mm4["m1"])
    diff_snp_sites(mm1["m1"], mm4["m1"])
    diff_snp_sites(mm1["m1"], mm4["m2"])
    diff_snp_sites(mm1["m2"], mm4["m1"])
    diff_snp_sites(mm1["m2"], mm4["m2"])


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


def common_pos_haplotypes(m1, m4, same_pos):
    h1 = m1.get_haplotypes()
    h4 = m4.get_haplotypes()
    new_h1 = {"m1": {}, "m2": {}}
    new_h4 = {"m1": {}, "m2": {}}
    #same_pos = find_common_poses(m1, m4)
    for pos in same_pos:
        new_h1["m1"][pos] = h1["m1"][pos]
        new_h1["m2"][pos] = h1["m2"][pos]
        new_h4["m1"][pos] = h4["m1"][pos]
        new_h4["m2"][pos] = h4["m2"][pos]
    return new_h1, new_h4


if __name__ == "__main__":
    SNPS = load_objects("../data/chr19_snps.obj")
    READS = load_objects("../data/chr19_reads_matched.obj")
    #IR_READS = load_objects("../data/chr19_reads_ir.obj")
    #DUMMY_READS = load_objects("../data/dummy_reads.obj")
    #DUMMY_SNPS = load_objects("../data/dummy_snps.obj")
    #print(len(DUMMY_READS))
    #print(len(DUMMY_SNPS))

    # gene DNMT1, chr19: 10133345 - 10231286 bp

    m1 = run_model(SNPS, READS, 100)
    #mm1 = m1.cluster_reads(IR_READS)
    m4 = run_model(SNPS, READS, 100)
    #mm4 = m4.cluster_reads(IR_READS)
    #compare_clustering(m1, m4, IR_READS)
    #compare_models(m1, m4)
    h1 = m1.get_haplotypes()
    h4 = m4.get_haplotypes()
    new_h1 = {"m1": {}, "m2": {}}
    new_h4 = {"m1": {}, "m2": {}}
    same_pos = find_common_poses(m1, m4)
    for pos in same_pos:
        new_h1["m1"][pos] = h1["m1"][pos]
        new_h1["m2"][pos] = h1["m2"][pos]
        new_h4["m1"][pos] = h4["m1"][pos]
        new_h4["m2"][pos] = h4["m2"][pos]
    compare_models(new_h1, new_h4)
