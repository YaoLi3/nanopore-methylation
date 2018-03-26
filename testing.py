# ! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 28/02/2018
"""
from nanoporereads import *
from haplotypes import *
from rawsignal import *
from snps import *

#############
#  Testing  #
#############
if __name__ == "__main__":
    """Nanopore reads data"""
    extract_fastq("/dict/fast5_folder/")
    # Use Minimap2 map reads to reference, get SAM file
    # DATA = NanoporeReads("data/chr19_merged.sam", "19")
    # DATA.get_reads()  # 45946 reads
    # overlap = DATA.find_imprinted(ImprintedRegions("data/ip_gene_pos.txt").get_regions(), 0, True, "data/find_imprinted_result.txt")
    #  375 reads

    """reads & snps data"""
    snps_data = load_VCF("data/chr19.vcf")  # list, 50746 SNPs on chr19
    reads_data = process_all_reads("data/find_imprinted_result.txt", snps_data)  # list, 302 reads

    # snps = []
    # for read in reads_data:
    # for snp in read.snps:
    # if snp not in snps:
    # snps.append(snp)

    """train HMM"""


    def read_log_likelihood(t, state, r, SNPs):
        """
        Calculate the log value of probability of a read observed in a model.
        :param t: theta, emission probs, 1-D numpy array
        :param state: one of the model, int 0 or 1.
        :param r: read
        :return log value of likelihood of the read in given model
        """
        read_llhd = 0
        observations = ["A", "G", "C", "T"]
        for snp in r.snps:
            base_llhd = t[SNPs.index(snp), state, observations.index(r.get_base(snp.pos))]
            read_llhd += math.log2(base_llhd)
        return read_llhd


    def assign_reads(emission_prob, all_reads, snps):
        # reset assignments
        m0_reads = []
        m1_reads = []
        m0_llhd = []
        m1_llhd = []
        m0_loci = []
        m1_loci = []

        for read in all_reads:
            l0 = read_log_likelihood(emission_prob, 0, read, snps)
            l1 = read_log_likelihood(emission_prob, 1, read, snps)
            if l0 > l1:
                m0_reads.append(read)  # save read objects
                m0_llhd.append(l0)  # save read likelihood log value
                for snp in read.snps:
                    if snp not in m0_loci:
                        m0_loci.append(snp)
            elif l0 < l1:
                m1_reads.append(read)
                m1_llhd.append(l1)
                for snp in read.snps:
                    if snp not in m1_loci:
                        m1_loci.append(snp)
        return m0_reads, m1_reads, m0_loci, m1_loci


    def base_list(m0_reads):
        bases = []
        for read in m0_reads:
            base = read.get_bases()
            bases.append(base)


    def maximize_likelihood_for_each_snp(theta, state, m_reads, m_loci, SNPs):
        """
        For only 1 state
        Model 0
        """
        # M-step (first half)
        loci_llhd = {}
        for locus in m_loci:
            theta_locus = theta[SNPs.index(locus), state,]  # a atcg distribution on this position

            locus_llhd = one_locus_llhd(theta_locus, locus, m_reads)
            loci_llhd[locus] = locus_llhd
            best_locus_atgc_distri = minimize(one_locus_llhd, theta_locus, args=(locus, m_reads), method="BFGS").x
            # print("best atgc emission prob on locus {} is {}".format(locus.pos, best_locus_atgc_distri))
            # print(best_locus_atgc_distri, theta_locus) # changes did happen on at least one base
            new_emission_matrix = update_emission_matrix(locus, SNPs, state, best_locus_atgc_distri, theta)
        return new_emission_matrix


    def update_emission_matrix(locus, SNPs, state, best_emission_on_one_locus, previous_emission_matirx):
        """
        For only one state.
        Model 0.
        """
        locus_in_matrix = SNPs.index(locus)
        emission_matrix = copy.deepcopy(previous_emission_matirx)
        emission_matrix[locus_in_matrix, state] = best_emission_on_one_locus
        return emission_matrix


    def one_locus_llhd(locus_theta, locussnp, m0_reads):
        reads_on_locus = []
        observation = ["A", "G", "C", "T"]
        locus_llhd = 0
        for read in m0_reads:
            if read.start <= locussnp.pos <= read.end:
                reads_on_locus.append(read)
        try:
            for read in reads_on_locus:
                base_on_locus = read.get_base(locussnp.pos)
                base_prob = math.log2(locus_theta[observation.index(base_on_locus)])
                # print("on {} locus, the likelihood of read {} has base {} is {}"
                # .format(locussnp.pos, read.id, base_on_locus, base_prob))
                locus_llhd += -base_prob
        except ValueError:
            pass
        return locus_llhd


    def if_sum_one(locus_theta):
        if sum(locus_theta) == 1:
            return True
        else:
            raise ValueError("Probabilities on this position do not sum to 1.")


    def if_array_sum_one(theta):
        n = 0;
        m = 0
        for (x, y, z), base_prob in np.ndenumerate(theta):
            if sum(theta[x, y,]) != 1:  # index is right
                n += 1
                m += 1
                print(x, y)
            else:
                n += 1
                pass
        print("final", n, m)


    def process(matrix, m0, m0_pos, m1, m1_pos, snps_data):
        matrix_0 = maximize_likelihood_for_each_snp(matrix, 0, m0, m0_pos, snps_data)
        matrix_1 = maximize_likelihood_for_each_snp(matrix, 1, m1, m1_pos, snps_data)
        m0_new, m1_new, m0_pos_new, m1_pos_new = assign_reads(matrix_1, reads_data, snps_data)
        return matrix_1, m0_new, m1_new, m0_pos_new, m1_pos_new


    h = Hmm(snps_data, reads_data, ["A", "G", "C", "T"], ["P", "M"])

    try:
        matrix_0 = h.init_emission_matrix()
        m0, m1, m0_pos, m1_pos = assign_reads(matrix_0, reads_data, snps_data)

        matrix_0_1 = maximize_likelihood_for_each_snp(matrix_0, 0, m0, m0_pos, snps_data)
        matrix_1 = maximize_likelihood_for_each_snp(matrix_0, 1, m0, m0_pos, snps_data)
        m01, m11, m0_pos1, m1_pos1 = assign_reads(matrix_1, reads_data, snps_data)

        matrix_0_2 = maximize_likelihood_for_each_snp(matrix_1, 0, m01, m0_pos1, snps_data)
        matrix_2 = maximize_likelihood_for_each_snp(matrix_1, 1, m01, m0_pos1, snps_data)
        m02, m12, m0_pos2, m1_pos2 = assign_reads(matrix_2, reads_data, snps_data)

        matrix_0_3 = maximize_likelihood_for_each_snp(matrix_2, 0, m01, m0_pos2, snps_data)
        matrix_3 = maximize_likelihood_for_each_snp(matrix_2, 1, m01, m0_pos2, snps_data)
        m03, m13, m0_pos3, m1_pos3 = assign_reads(matrix_3, reads_data, snps_data)

        matrix_0_4 = maximize_likelihood_for_each_snp(matrix_3, 0, m01, m0_pos3, snps_data)
        matrix_4 = maximize_likelihood_for_each_snp(matrix_3, 1, m01, m0_pos3, snps_data)
        m04, m14, m0_pos4, m1_pos4 = assign_reads(matrix_4, reads_data, snps_data)

        matrix_5, m05, m15, m0_pos5, m1_pos5 = process(matrix_4, m04, m0_pos4, m14, m1_pos4, snps_data)
        print(len(m0), len(m1))
        print(len(m01), len(m11))
        print(len(m02), len(m12))
        print(len(m03), len(m13))
        print(len(m04), len(m14))
        print(len(m05), len(m15))

        A = []; B = []; T = []
        for read in m05:
            a = 0; b = 0
            for snp in read.snps:
                if snp.gt == "0|1":
                    a += 1
                elif snp.gt == "1|0":
                    b += 1
            print(a, b)
            if a > b:
                A.append(read)
            elif a < b:
                B.append(read)
            else:
                T.append(read)

        A1 = [];
        B1 = [];
        T1 = []
        for read in m15:
            a = 0;
            b = 0
            for snp in read.snps:
                if snp.gt == "0|1":
                    a += 1
                elif snp.gt == "1|0":
                    b += 1
            if a > b:
                A1.append(read)
            elif a < b:
                B1.append(read)
            else:
                T1.append(read)

        print("Model 0: A{} B{} T{}".format(len(A), len(B), len(T)))
        print("Model 1: A{} B{} T{}".format(len(A1), len(B1), len(T1)))
        """
        a = 0; b = 0; c = 0; d = 0
        ref0 = 0; ref1 = 0; alt0 = 0; alt1 = 0; other0 = 0; other1 = 0
        r0 = []; a0 = []; r1 = []; a1 = []
        for index, read in enumerate(m05):
            for snp in read.snps:
                base = read.get_base(snp.pos)
                if base == snp.ref:
                    ref0 += 1
                    r0.append(snp)
                elif base == snp.alt:
                    alt0 += 1
                    a0.append(snp)
                else:
                    other0 += 1
        for snp in r0:
            if snp.gt == "0|1":
                a += 1
            elif snp.gt == "1|0":
                b += 1
        for snp in a0:
            if snp.gt == "0|1":
                b += 1
            elif snp.gt == "1|0":
                a += 1

        for read in m15:
            for snp in read.snps:
                base = read.get_base(snp.pos)
                if base == snp.ref:
                    ref1 += 1
                    r1.append(snp)
                elif base == snp.alt:
                    alt1 += 1
                    a1.append(snp)
                else:
                    other1 += 1
        for snp in r1:
            if snp.gt == "0|1":
                c += 1
            elif snp.gt == "1|0":
                d += 1
        for snp in a1:
            if snp.gt == "0|1":
                d += 1
            elif snp.gt == "1|0":
                c += 1

        print("Model0: A{}, B{}".format(a, b))
        print("Model1: A{}, B{}".format(c, d))
        """
        #n = 0
        #while n <= 100:
            #matrix_first, m0_first, m1_first, m0_pos_f, m1_pos_f = process(matrix_5, m05, m0_pos5, m15, m1_pos5,
                                                                           #snps_data)
            #matrix, m0, m1, m0_pos, m1_pos = process(matrix_first, m0_first, m0_pos_f, m1_first, m1_pos_f, snps_data)
            #print(len(m0), len(m1))
            #n += 1
    except ValueError:
        pass


    def one_init_a_time():
        init_theta = h.init_emission_matrix()
        m0, m1, m0_pos, m1_pos = assign_reads(init_theta, reads_data, snps_data)
        matrix, m0, m1, m0_pos, m1_pos = process(matrix_first, m0_first, m0_pos_f, m1_first, m1_pos_f, snps_data)



    """raw signals"""
    # raws = get_raw_dirc("/shares/coin/yao.li/data/basecall_pass/", "/shares/coin/yao.li/signal/", overlap)
    # h1, h2 = find_haplotype(raws, haplotypes)
