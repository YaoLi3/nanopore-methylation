#! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 08/02/2018
"""
import numpy as np
import editdistance


#############################
#  Human Imprinted Regions  #
#############################
class ImprintedRegions:
    """
    Define the imprinted regions on human genome (EnSembl hg38).
    Save the data in a dictionary.
    """

    def __init__(self, filename):
        """
        :param filename: (string) path of the file stores the
                         information of certain imprinted regions.
        """
        self.regions = {}
        file = open(filename, "r")
        for line in file:
            line = line.strip().split("\t")
            self.regions[line[0]] = line[1:]
        file.close()

    def get_regions(self):
        """
        :return: (dictionary) key = gene names
        """
        return self.regions


##################
# Nanopore reads #
##################
class NR:
    def __init__(self, id, seq, chr, poses):
        self.id = id
        self.raw_signal = None
        self.sequence = seq
        self.chr = chr
        self.poses = poses
        self.overlap_region = None
        self.gene = ""
        self.overlap_fastq = ""
        self.overlap_pos = (0, 0)
        self.haplotype = ""
        self.snps = []


def load_Sam(samfile, chrom):
    data = []
    file = open(samfile, "r")
    for line in file:
        if line.startswith("@"):  # ignore headers
            pass
        else:
            line = line.split("\t")
            rname = line[2]
            if rname == chrom:  # some READs may map to other chromosomes
                start = int(line[3])
                seq = line[9]
                seq_len = len(seq)
                end = (start + seq_len)
                nr = NR(line[0], seq, rname, (start, end))
                data.append(nr)
    file.close()
    return data


class NanoporeReads:
    """
    Discover reads that are overlapped with human imprinted regions.
    From one sam file.
    """

    def __init__(self, samfile, chrom):
        """
        Open a sam file, extract data of reads
        mapped to the right chromosome on reference genome.
        :param samfile: (string) path of external sam file
        :param chrom: (string) chromosome number
        """
        self.samfile = samfile  # store filename as an attribute]
        try:
            if 1 <= int(chrom) <= 22:
                self.chrom = "chr" + chrom
        except ValueError:
            if chrom in ["XxYx"]:
                self.chrom = "chr" + chrom.lower()
            else:
                raise Exception("Invalid chromosome number.")

        file = open(self.samfile, "r")
        self.data = {}
        self.reads = {}
        self.overlap = {}
        for line in file:
            if line.startswith("@"):  # ignore headers
                pass
            else:
                line = line.split("\t")
                self.data[line[0]] = line[1:]
                # calculate the region of a READs mapped on to the reference genome
                start = int(line[3])
                seq = line[9]
                seq_len = len(seq)
                end = (start + seq_len)
                rname = line[2]
                if rname == self.chrom:  # some READs may map to other chromosomes
                    self.reads[line[0]] = (start, end, rname, seq)
        file.close()

    def get_data(self):
        """
        :return: (dictionary) all data from the sam file as a dictionary.
        """
        return self.data

    def get_reads(self):
        """
        :return: (dictionary) information of names,
                 start and end positions of reads.
        """
        return self.reads

    def get_sequence(self, read_name):
        """
        :param read_name: (string) Nanopore basecalling ID
        :return: (string) sequence of the given READs
        """
        return self.data[read_name][8]

    def search_reads(self, reads_names):
        """
        Search reads based on their given QNAMEs.
        :param reads_names: (string) Nanopore basecalling ID
        :return: boolean value
        """
        for qname in reads_names:
            try:
                self.data.get(qname)
            except ValueError:
                print("This READs does not exist.")
        return True

    def find_imprinted(self, regions, thrhld, save_file=False, file=None):
        """
        For a given sam file, find out if there is any READs in the file
        is located in human genetic imprinted regions.
        :param regions: (dictionary) positions of human imprinted regions on reference genome
        :param thrhld: (float) a number between 0 and 1. portion? does not matter
        :param save_file: (bool) if want to save the result to a txt file
        :param file: (string) path of the file to save results
        :return: (dictionary) key = imprinted reads ID
                            values: overlapped imprinted gene name
                            which end of the READs located in the imprinted region
                            positions of overlap segment refer to the original READs
                            positions of overlap segment refer to the reference genome
                            threshold been used
        """
        self.overlap = {}
        gene = {}; b = 0; s = 0; e = 0; n = 0; cnt = 0
        for j in regions:  # j = gene name
            c = 0
            start = int(regions[j][0])
            end = int(regions[j][1])
            r_range = range(start, end + 1)
            chrom = regions[j][2]  # chr is just a number(str)
            for i in self.reads:  # i = READs ID
                # make sam id and fast5 id the same format
                if i.find("_Basecall_1D_template"):
                    read_id = i.replace("_Basecall_1D_template", "")
                else:
                    read_id = i.replace("_Basecall_Alignment_template", "")

                pos1, pos2, rname, seq = self.reads[i]  # pos1 & pos2 are int
                if 0 <= thrhld <= 1:
                    min_coverage = thrhld * (pos2 - pos1)
                else:
                    print("t is a float between 0 and 1.")
                    return
                if chrom == rname.replace("chr", ""):
                    if pos1 in r_range and pos2 in r_range:
                        b += 1; cnt += 1
                        l = pos2 - pos1
                        # if both ends of the READs located in an imprinted region
                        if min_coverage <= l:
                            n += 1; c += 1
                            gene[j] = c
                            self.overlap[read_id] = [j, self.chrom, "both ends",
                                                     (1, l + 1),
                                                     (pos1, pos2),
                                                     thrhld, seq]
                    elif pos1 in r_range and pos2 not in r_range:
                        s += 1; cnt += 1
                        l1 = end - pos1  # overlapped length
                        if min_coverage <= l1:
                            n += 1; c += 1
                            gene[j] = c
                            self.overlap[read_id] = [j, self.chrom, "start pos",
                                                     (1, l1 + 1),
                                                     (end - l1, end),
                                                     thrhld, seq]
                    elif pos2 in r_range and pos1 not in r_range:
                        e += 1; cnt += 1
                        l2 = pos2 - start
                        if min_coverage <= l2:
                            n += 1; c += 1
                            gene[j] = c
                            self.overlap[read_id] = [j, self.chrom, "end pos",
                                                     (pos2 - pos1 - l2, pos2 - pos1),
                                                     (start, start + l2), thrhld, seq]
        # Save results into a txt file
        if save_file:
            file = open(file, "w")
            file.write(
                "Read_ID\t\tImprinted_Gene\t\tChromosome\t\tInfo\t\tPos_On_Read"
                "\t\tPos_On_Ref_Genome\t\tIR_Length_Threshold\n\n")
            for id in self.overlap:
                file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n\n".format(id, self.overlap[id][0],
                                                                       self.overlap[id][1],
                                                                       self.overlap[id][2],
                                                                       self.overlap[id][3],
                                                                       self.overlap[id][4],
                                                                       self.overlap[id][5],
                                                                       self.overlap[id][6]))
            file.write("\n{} reads pass the threshold.\n".format(n))
            for name in gene:
                file.write("{} in gene {}\n".format(gene[name], name))
            file.write(
                "\nTotal {} reads have both ends located in "
                "an imprinted region.\n{} reads start position"
                " mapped to imprinted region.\n{} reads end position"
                " mapped to imprinted region.\n".format(b, s, e))
            file.write("\nTotal {} chr19 reads have overlap "
                       "with chr19 imprinted regions.\n".format(cnt))
            file.close()

    def get_imprinted_reads(self):
        """
        :return: (list) imprinted Nanopore reads IDs
        """
        return self.overlap

    def search_imprinted_read(self, ID):
        return self.overlap[ID]

    def get_matrix(self):
        """
        Calculate the minimum edit distance between two reads.
        :return: (numpy narray) pairwise distance matrix
        """
        if not self.overlap == []:
            dist_matrix = np.zeros((len(self.overlap), len(self.overlap)))
            a = 0
            for i in self.overlap:
                seq1 = self.get_sequence(i)
                b = 0
                for j in self.overlap:
                    seq2 = self.get_sequence(j)
                    dist_matrix[a, b] = editdistance.eval(seq1, seq2)
                    b += 1
                a += 1
            return dist_matrix
        else:
            print("The list of imprinted reads is empty.")

    def get_read_len(self, read):
        """for id in self.reads:
                read_len = self.get_read_len(self.reads[key])"""
        return len(read[3])


###################
# Human SNPs data #
###################
class SNPs:
    """
    Human SNPs data.
    chrom, pos, id, ref, alt, qual, filter, info, gt
    (ignore all indel variants)
    How many of them locate in nanopore reads? in human imprinted regions?

    e.g.
    chr19 294525 . A C 0 PASS KM=8;KFP=0;KFF=0;MTD=bwa_freebayes,bwa_gatk,bwa_platypus,isaac_strelka GT 1|0
    """

    def __init__(self, chr, id, pos, ref, alt, gt):
        self.chrom = chr
        self.id = id
        self.pos = int(pos)
        self.ref = ref
        self.alt = alt
        self.mut = ref + alt
        self.gt = gt
        self.type = ""
        self.detect_type()

    def detect_type(self):
        """
        Determine a SNP type.
        """
        TRANSITIONS = ["AG", "CT"]
        TRANSVERSIONS = ["AC", "AT", "CG", "GT"]

        if len(self.ref) > 1 or len(self.alt) > 1:
            self.type = "indel"

        elif (self.ref in TRANSITIONS[0] and self.alt in TRANSITIONS[0]) \
                or (self.ref in TRANSITIONS[1] and self.alt in TRANSITIONS[1]):
            self.type = "transition"

        else:
            for conbo in TRANSVERSIONS:
                if self.ref != self.alt and self.ref in conbo and self.alt in conbo:
                    self.type = "transversion"

    def __str__(self):
        return "{}: {}\tREF:{}, ALT:{}\tTYPE:{}".format(self.chrom, self.pos, self.ref, self.alt, self.type)

    def __eq__(self, other):
        """Override the default Equals behavior"""
        return self.chrom == other.chrom and self.pos == other.pos and self.mut == other.mut

    def __ne__(self, other):
        """Override the default Unequal behavior"""
        return self.chrom != other.chrom or self.pos != other.pos or self.mut != other.mut


def load_VCF(vcf_file, count_index=False):
    """
    Read a VCF file and return he data in it.
    :param vcf_file: (string) VCF file name
    :return: all_snps: (list) SNPs instances
    """
    try:
        f = open(vcf_file, "r")
        all_snps = []
        for line in f:
            if line.startswith("chr"):
                chrom, pos, id, ref, alt, qual, \
                filter, info, format, gt \
                    = line.strip().split("\t")
                a, b = gt.split("|")
                if chrom == "chr19" and a != b:  # only use heter snps
                    snp = SNPs(chrom, id, pos, ref, alt, gt)
                    if not snp.type == "indel":
                        all_snps.append(snp)
        f.close()
        return all_snps
    except ValueError:
        raise RuntimeError("Not the right values to unpack.")
    except IOError:
        raise IOError("This vcf file is not available.")


def count_snps(all_snps):
    """
    :param all_snps:
    :return:
    """
    MUTATIONS = {"AG": 0, "AC": 0, "AT": 0, "CA": 0, "CT": 0, "CG": 0,
                 "TA": 0, "TC": 0, "TG": 0, "GA": 0, "GC": 0, "GT": 0}
    REF = {"A": 0, "T": 0, "G": 0, "C": 0}
    ALT = {"A": 0, "T": 0, "G": 0, "C": 0}

    for snp in all_snps:
        for mutation in MUTATIONS:
            if snp.mut == mutation:
                MUTATIONS[mutation] += 1

    for snp in all_snps:
        for r in REF:
            if snp.ref == r:
                REF[r] += 1
        for a in ALT:
            if snp.alt == a:
                ALT[a] += 1
    return MUTATIONS, REF, ALT


def get_snp_prob(all_snps, snp):
    """
    Given a SNPs, what's the frequency of happening
    its mutation (ref -> alt)? (ignore positions?)
    :param all_snps: list of SNPs objects
    :param snp: (SNPs) an instance
    :return: probability (float) calculated based on give vcf data.
    """
    counts, ref, alt = count_snps(all_snps)
    total = 0
    ref_t = 0
    alt_t = 0
    p = 0
    p_ref = 0
    p_alt = 0

    for key in counts:
        total += counts[key]
    if snp.mut in counts:
        p = counts[snp.mut] / total
    else:
        print("Wrong ref or alt value.")

    for r in ref:
        ref_t += ref[r]
    if snp.ref in ref:
        p_ref = ref[snp.ref] / ref_t
    else:
        print("Wrong ref or alt value.")

    for a in alt:
        alt_t += alt[a]
    if snp.alt in alt:
        p_alt = alt[snp.alt] / alt_t
    else:
        print("Wrong ref or alt value.")

    return p_ref * p, p_alt * p


def count_types(all_snps):
    """
    Count numbers of each type of SNPs.
    :return: lists
    """
    indels = []
    transitions = []
    transversions = []
    for snp in all_snps:
        if snp.type == "indel":
            indels.append(snp)
        elif snp.type == "transition":
            transitions.append(snp)
        else:
            transversions.append(snp)
    return indels, transitions, transversions


def get_positions(all_snps):
    pos = []
    for snp in all_snps:
        pos.append(snp.pos)
    return pos


########################
# Target READs segments #
########################
class OverlappedRead:
    """
    A READs object has 5 attributes:
    its sequencer id,
    chromosome number,
    start position on reference genome (gh38),
    end position on reference genome,
    and its SNPs objects.

    In this step, we only deal with nanopore reads that have at least one base
    overlapped with the known human imprinted gene regions.

    And to find SNPs located inside each reads using known human SNPs data.
    """

    def __init__(self, id, chr, pos1, pos2, seq):
        self.id = id
        self.chrom = chr
        self.start = int(pos1)
        self.end = int(pos2)
        self.seq = seq
        self.snps = []
        self.base = []
        self.state = ""

    def detect_snps(self, SNPs_data):
        """
        Find snps in one READs.
        :param SNPs_data: (list) class SNPs objects
        """
        for snp in SNPs_data:  # each SNPs instance, attr: chrom, id, pos, ref, alt, gt
            if snp.chrom == self.chrom and self.start <= snp.pos <= self.end:  # if position
                self.snps.append(snp)

    def get_base(self, pos):
        try:
            index = int(pos) - int(self.start)
            return self.seq[index]
        except IndexError:
            print("snp pos not right.")

    def get_bases(self):
        # TODO: make sure the bases are right. positions!
        """
        Get bases on SNP positions on READs sequence.
        :return: list of string, bases, ATCG
        """
        for snp in self.snps:
            index = int(snp.pos) - int(self.start)
            self.base.append(self.seq[index])
        return self.base

    def get_read_data(self):
        return self.id, self.snps, self.base

    def set_state(self, state):
        """Set hidden markov STATEs for a READs."""
        self.state = state

    def __str__(self):
        return "{}: {},{} \tSNPs:{}".format(self.chrom, self.start, self.end, len(self.snps))

    def __eq__(self, other):
        """Override the default Equals behavior"""
        return self.id == other.id and self.state == other.state

    def __ne__(self, other):
        """Override the default Unequal behavior"""
        return self.id != other.id or self.snps != other.snps or self.state != other.state


def process_all_reads(read_file, SNPs_data):
    """
    1.Nanopore reads data.
    2.SNPs data.

    :param read_file:
    :param SNPs_data:
    :return: (list) : OverlappedReads instances (total 375 of them)
    """
    f = open(read_file, "r")
    all_reads = []
    for line in f:
        line = line.strip().split("\t")
        if len(line) == 8:
            id, gene, chr, info, read_pos, ref_pos, thrhld, seq = line
            ref_pos = ref_pos.strip("()").split(",")
            read = OverlappedRead(id, chr, ref_pos[0], ref_pos[1][1:], seq)
            read.detect_snps(SNPs_data)
            if not read.snps == []:
                read.get_bases()
                all_reads.append(read)
    f.close()
    return all_reads


########################################################
# Manually find variants between ref and nanopore seqs #
########################################################
def read_ref_genome(fn, chr):
    """
    :param fn:
    :param chr: (int)
    :return:
    """
    chrs = ["chr1, chr19"]
    l = open(fn).readlines()
    region = []
    for index, line in enumerate(l):
        # line = line.strip().split()
        # print(line)
        if line.startswith(">chr{} ".format(chr)):
            region.append(index + 1)
        if line.startswith(">chr{} ".format(chr + 1)):  # need to change this part
            region.append(index)
    seq = ("".join(l[region[0]: region[1]])).replace("\n", "")  # 58617616 base pairs, chr19
    return seq


def compare_seqs(nano_seq, ref_seq):
    """
    :param nano_seq:
    :param ref_seq:
    :return:(dict) SNPs/same: key = position on the reads, value = (ref, alt)
    """
    snp = {}
    same = {}
    if len(nano_seq) != len(ref_seq):
        raise ValueError("Sequences don't have the same length.")
    else:
        for pos in range(len(nano_seq)):
            if nano_seq[pos] != ref_seq[pos]:
                snp[pos] = (nano_seq[pos], ref_seq[pos])
            elif nano_seq[pos] == ref_seq[pos]:
                same[pos] = (nano_seq[pos], ref_seq[pos])
        return snp, same


def find_snps_in_aln(ip_reads, ref_genome, chrom="chr19"):
    """
    :param chrom:
    :param ip_reads: (dict) imprinted nanopore reads results, id is key, ref positions is [3]
    :param ref_genome:
    :return:
    """
    nanopore_snp = {}
    read = {}
    for id in ip_reads:  # iterate through each nanopore READs sequence
        # get where the nanopore READs mapped to the reference genome
        ref_start, ref_end = (int(ip_reads[id][3][0]),
                              int(ip_reads[id][3][1][1:]))
        read_start, read_end = (int(ip_reads[id][2][0]),
                                int(ip_reads[id][2][1][1:]))

        aln_type = {}
        # get nanopore sequence and its correlated reference genome sequence to compare
        nanopore_seq = ip_reads[id][5][read_start - 1: read_end - 1]
        ref_seq = ref_genome[ref_start - 1: ref_end - 1]

        snp, same = compare_seqs(nanopore_seq, ref_seq)

        aln_type["SNPs"] = snp
        aln_type["same"] = same
        read[id] = aln_type
    nanopore_snp[chrom] = read
    return nanopore_snp