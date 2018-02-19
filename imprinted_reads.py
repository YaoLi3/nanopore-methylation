"""
__author__ = Yao LI
__email__ = liyaoo1012@163.com
__date__ = 08/02/2018
"""
import h5py
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

    def getRegions(self):
        """
        :return: (dictionary) key = gene names
        """
        return self.regions


####################
#  Nanopore reads  #
####################
class NanoporeReads:
    """
    Discover reads that are overlapped with human imprinted regions.
    From one sam file.
    """

    def __init__(self, samfile, chrom="chr"):
        """
        Open a sam file, extract data of reads
        mapped to the right chromosome on reference genome.
        :param samfile: (string) path of external sam file
        :param chrom: (string) chromosome name(s)
        """
        self.samfile = samfile  # store filename as an attribute
        self.chrom = chrom
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
                # calculate the region of a read mapped on to the reference genome
                start = int(line[3])
                seq_len = len(line[9])
                end = (start + seq_len)
                rname = line[2]
                #if rname = self.chrom
                if rname == "chr19":  # only extract reads mapped to chr19 (ref)
                    self.reads[line[0]] = (start, end, rname)
        file.close()

    def getData(self):
        """
        :return: (dictionary) all data from the sam file as a dictionary.
        """
        return self.data

    def getReads(self):
        """
        :return: (dictionary) information of names,
                 start and end positions of reads.
        """
        return self.reads

    def getSequence(self, read_name):
        """
        :param read_name: (string) Nanopore basecalling ID
        :return: (string) sequence of the given read
        """
        return self.data[read_name][8]

    def searchReads(self, reads_names):
        """
        Search reads based on their given QNAMEs.
        :param reads_names: (string) Nanopore basecalling ID
        :return: boolean value
        """
        for qname in reads_names:
            try:
                self.data.get(qname)
            except ValueError:
                print("This read does not exist.")
        return True

    def findImprinted(self, regions, thrhld, save_file=False, file=None):
        """
        For a given sam file, find out if there is any read in the file
        is located in human genetic imprinted regions.
        :param regions: (dictionary) positions of human imprinted regions on reference genome
        :param thrhld: (float) a number between 0 and 1. portion? does not matter
        :param save_file: (bool) if want to save the result to a txt file
        :param file: (string) path of the file to save results
        :return: (dictionary) key = imprinted reads ID
                            values: overlapped imprinted gene name
                            which end of the read located in the imprinted region
                            positions of overlap segment refer to the original read
                            positions of overlap segment refer to the reference genome
                            threshold been used
        """
        self.overlap = {}
        gene = {}
        b = 0
        s = 0
        e = 0
        n = 0
        cnt = 0
        for j in regions:  # j = gene name
            c = 0
            start = int(regions[j][0])
            end = int(regions[j][1])
            r_range = range(start, end + 1)
            chrom = regions[j][2]  # chr is just a (str)number
            for i in self.reads:  # i = read ID
                pos1, pos2, rname = self.reads[i]  # pos1 & pos2 are int
                if 0 <= thrhld <= 1:
                    min_coverage = thrhld * (pos2 - pos1)
                else:
                    print("t is a float between 0 and 1.")
                    return
                if chrom == rname.replace("chr", ""):
                    if pos1 in r_range and pos2 in r_range:
                        b += 1
                        cnt += 1
                        l = pos2 - pos1
                        # if both ends of the read located in an imprinted region
                        if min_coverage <= l:
                            n += 1
                            c += 1
                            gene[j] = c
                            self.overlap[i] = [j, "both ends",
                                               (1, l + 1),
                                               (pos1, pos2),
                                               thrhld]
                    elif pos1 in r_range and pos2 not in r_range:
                        s += 1
                        cnt += 1
                        l1 = end - pos1  # overlapped length
                        if min_coverage <= l1:
                            n += 1
                            c += 1
                            gene[j] = c
                            self.overlap[i] = [j, "start pos",
                                               (1, l1 + 1),
                                               (end - l1, end),
                                               thrhld]
                    elif pos2 in r_range and pos1 not in r_range:
                        e += 1
                        cnt += 1
                        l2 = pos2 - start
                        if min_coverage <= l2:
                            n += 1
                            c += 1
                            gene[j] = c
                            self.overlap[i] = [j, "end pos",
                                               (pos2 - pos1 - l2, pos2 - pos1),
                                               (start, start + l2), thrhld]
        # Save to a txt file
        if save_file:
            file = open(file, "w")
            file.write(
                "Read_ID\t\tImprinted_Gene\t\tInfo\t\tPos_On_Read"
                "\t\tPos_On_Ref_Genome\t\tIR_Length_Threshold\n\n")
            for id in self.overlap:
                file.write("{}\t{}\t{}\t{}\t{}\t{}\n\n".format(id, self.overlap[id][0],
                                                               self.overlap[id][1],
                                                               self.overlap[id][2],
                                                               self.overlap[id][3],
                                                               self.overlap[id][4]))
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

        return self.overlap

    def getImprintedReads(self):
        """
        :return: (list) imprinted Nanopore reads IDs
        """
        return self.overlap

    def searchIR(self, ID):
        return self.overlap[ID]

    def getMatrix(self):
        """
        Calculate the minimum edit distance between two reads.
        :return: (numpy narray) pairwise distance matrix
        """
        if not self.overlap == []:
            dist_matrix = np.zeros((len(self.overlap), len(self.overlap)))
            a = 0
            for i in self.overlap:
                seq1 = self.getSequence(i)
                b = 0
                for j in self.overlap:
                    seq2 = self.getSequence(j)
                    dist_matrix[a, b] = editdistance.eval(seq1, seq2)
                    b += 1
                a += 1
            return dist_matrix
        else:
            print("The list of imprinted reads is empty.")


def extract_fastq(name):
    """
    Extract fastq sequence from a fast5 file.
    :param name: fast5 file name
    :return: (list) a fastq sequence
    """
    f = h5py.File(name, "r")
    seq = f['Analyses']['Basecall_1D_001']['BaseCalled_template']['Fastq'].value
    if seq != "":
        fastq_file = open("%s.fastq" % name, "wb")
        fastq_file.write(seq)
        fastq_file.close()

    seq = []
    f = open("%s.fastq" % name, "r")
    for line in f:
        seq.append(line)
    f.close()
    return seq

