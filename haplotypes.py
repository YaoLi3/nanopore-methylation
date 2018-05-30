import numpy as np
import math
import collections  # to sort dictionaries


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
            try:
                t = self.get_emission()
                base_llhd = t[snp.id, state, self.observations.index(read.bases[snp.id])]
                read_llhd += math.log(base_llhd)  # math domain error, base_llhd = -8.84548072603065
            except ValueError:
                continue
        return read_llhd

    def assign_reads(self, reads):
        """Calculate the parental-maternal likelihood of the reads.
        Expectation step in EM algorithm.
        
        Args:
            reads: list of reads.
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
                    self.alleles["m1"][self.SNPs[snp_id].pos] = read.bases[snp_id]
            elif pm_llhd[idx, 0] < pm_llhd[idx, 1]:
                self.read_assignments["m2"].append(idx)
                for snp_id in read.snps_id:
                    self.alleles["m2"][self.SNPs[snp_id].pos] = read.bases[snp_id]
        return pm_llhd

    @staticmethod
    def random_choose(l):
        """must be one-d, a list or an integer works. dictionaries not."""
        s = len(l) // 2
        return np.random.choice(l, size=s, replace=False)

    @staticmethod
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

    def update(self, reads, sr_dict, pm_llhd, hard=False, pseudo_base=1e-1):
        """Update the probability matrix
        
        Args:
            reads (List): List of reads.
            sr_dict(Dictionary): Dictionary of snp snp map.
            pm_llhd (Float): Matrix of shape [Read_num,2].
            hard (Boolean): If update in hard manner.
            pseudo_base(Float): Probability normalized by adding this pseudo_base, 
                e.diffs. p[_] = (count[_]+pseudo_base)/(count[_]+4*pseudo_base).
            restrict (Boolean): if selectively use read snps

        hard: count[self.observations.index(read.get_base(snp.pos)), pm_type] += 1
        fuzzy: count[self.observations.index(read.get_base(snp.pos)), 0] += (pm_llhd[read_id, 0] * 1)
        """
        for snp_id in sr_dict.keys():
            #snp = self.SNPs[snp_id]
            count = np.zeros((len(self.observations), 2))
            count[:] = pseudo_base
            for read_id in sr_dict[snp_id]:
                try:
                    read = reads[read_id]
                    if hard:
                        pm_type = np.argmax(pm_llhd[read_id, :])
                        count[self.observations.index(read.bases[snp_id]), pm_type] += 1  # snp_id or snp.id, same
                    else:
                        count[self.observations.index(read.bases[snp_id]), 0] += (pm_llhd[read_id, 0] * 1)
                        count[self.observations.index(read.bases[snp_id]), 1] += (pm_llhd[read_id, 1] * 1)
                except ValueError:  # TODO: no value error should occur
                    continue
                for state in range(len(self.STATEs)):
                    self.emission_probs[snp_id, state, :] = (count[:, state] / np.sum(count[:, state]))  # normalize

    def alter_update(self, reads, sr_dict, pm_llhd, pseudo_base=1e-1):
        """Only consider a part of SNPs site in a read.
        soft update manner"""
        snps = self.random_choose(list(sr_dict.keys()))
        for snp_id in snps:
            snp = self.SNPs[snp_id]
            count = np.zeros((len(self.observations), 2))
            count[:] = pseudo_base
            for read_id in sr_dict[snp_id]:
                try:
                    read = reads[read_id]
                    count[self.observations.index(read.bases[snp.id]), 0] += (pm_llhd[read_id, 0] * 1)  # add counts
                    count[self.observations.index(read.bases[snp.id]), 1] += (pm_llhd[read_id, 1] * 1)
                except ValueError:
                    continue
            for state in range(len(self.STATEs)):
                self.emission_probs[snp_id, state, :] = (count[:, state] / np.sum(count[:, state]))  # normalize

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
        """Return two lists of SNPs loci."""
        p1 = sorted(list(self.alleles["m1"].keys()))
        p2 = sorted(list(self.alleles["m2"].keys()))
        diff1 = []
        diff2 = []
        if len(p1) != len(p2):
            print("The number of SNPs in two states are different.")
            for snp in p1:
                if snp not in p2:
                    diff1.append(snp)
            print("These SNP loci only appear in state 1: {}".format(diff1))
            for snp in p2:
                if snp not in p1:
                    diff2.append(snp)
            print("These SNP loci only appear in state 2: {}".format(diff2))
        elif p1 == p2:
            print("SNPs loci generated by two states are the same.")
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


def run_model(snps, reads, iter_num, hard=False, updateAll=True):
    """Set iteration times, run a model.
    Alternative: iterate until theta converge."""
    model = HMM(snps, ["Parental", "Maternal"])
    model.init_emission_matrix()
    s = np.zeros((iter_num, 2))
    for _ in range(iter_num):
        sr_dict = model.snp_read_dict(reads)
        pm_llhd = model.assign_reads(reads)
        #print(model.emission_probs[9519,:])
        #print((np.sum(pm_llhd[:,1]),np.sum(pm_llhd[:,0])))
        s[_, 0] = np.sum(pm_llhd[:, 0])
        s[_, 1] = np.sum(pm_llhd[:, 1])
        if updateAll:
            model.update(reads, sr_dict, pm_llhd, hard=hard)
        else:
            model.alter_update(reads, sr_dict, pm_llhd)
    return model


def models_iterations(iter_times, snps, reads, hard=False):
    """Generate a markov model multiple times."""
    model = run_model(snps, reads, 10, hard=hard)
    for _ in range(iter_times):
        last_model = model
        model = run_model(snps, reads, 10, hard=hard)
        compare_models(last_model, model)


def compare_haplotypes(h11, h12):
    """
    Compare two haplotype dictionaries.
    :param h11: (dict) {snp_pos: base, snp_base: base}
    :param h12: (dict) {snp_pos: base, snp_base: base}
    :return: diff (dict) Missing, Un-match SNP sites and read bases.
    """
    diff = {"h1": {}, "h2": {}, "diff": {}}
    for key in h11:
        if key not in h12:
            diff["h1"][key] = h11[key]
        elif h11[key] != h12[key]:
            diff["diff"][key] = (h11[key], h12[key])
    for k in h12:
        if k not in h11:
            diff["h2"][k] = h12[k]
    return diff


def compare_models(m1, m2):
    """
    Compare haplotype clusters generate by two hidden markov models.
    :param m1: (HMM object) a HMM hidden markov model
    :param m2: (HMM) a hidden markov model
    :return: (None) print out changes between clusters. (addition, missing, different bases on the same SNP sites)
    """
    h11 = m1.get_alleles()["m1"]
    h12 = m1.get_alleles()["m2"]
    h21 = m2.get_alleles()["m1"]
    h22 = m2.get_alleles()["m2"]

    # h1 vs h1, h2 vs h2
    d1 = compare_haplotypes(h11, h21)
    d2 = compare_haplotypes(h12, h22)
    # h1 vs h2, h2 vs h1
    d3 = compare_haplotypes(h11, h22)
    d4 = compare_haplotypes(h12, h21)

    print("model1 cluster1 has {} alleles, model2 cluster1 has {} alleles.\n{} bases are different, "
          "{} bases only in model1 cluster1, {} bases only in model2 cluster1.".format(len(h11), len(h21),
                                                                                   len(d1["diff"]),
                                                                                   len(d1["h1"]),
                                                                                   len(d1["h2"])))
    print("model1 cluster2 has {} alleles, model2 cluster2 has {} alleles.\n{} bases are different, "
          "{} bases only in model1 state2, {} bases only in model2 state2.".format(len(h12), len(h22),
                                                                                   len(d2["diff"]),
                                                                                   len(d2["h1"]),
                                                                                   len(d2["h2"])))
    print("model1 cluster1 has {} alleles, model2 cluster2 has {} alleles.\n{} bases are different, "
          "{} bases only in model1 cluster1, {} bases only in model2 cluster2.".format(len(h11), len(h22),
                                                                                   len(d3["diff"]),
                                                                                   len(d3["h1"]),
                                                                                   len(d3["h2"])))
    print("model1 cluster2 has {} alleles, model2 cluster1 has {} alleles.\n{} bases are different, "
          "{} bases only in model1 cluster2, {} bases only in model2 cluster1.\n".format(len(h12), len(h21),
                                                                                   len(d4["diff"]),
                                                                                   len(d4["h1"]),
                                                                                   len(d4["h2"])))


def gold_standard(g_dict, m_dict):
    """compare model dict with this gold standard dict
    only use SNPs loci that been covered by Nanopore read in the model
    no need to use all the SNPs from the VCF file"""
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
