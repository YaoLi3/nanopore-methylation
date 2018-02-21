"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 20/02/2018
"""
import os
import h5py


def extractFastq(fast5_fn):
    """
    Extract fastq sequence from a given fast5 file,
    and save the fastq sequence in a file with the same name.
    :param fast5_fn: (string) fast5 file name
    """
    f = h5py.File(fast5_fn, "r")
    seq = f['Analyses']['Basecall_1D_001']['BaseCalled_template']['Fastq'].value
    if seq != "":
        fastq_file = open("%s.fastq"% fast5_fn.replace(".fast5", ""), "wb")
        fastq_file.write(seq)
        fastq_file.close()


def readFastq(fastq_fn):
    """
    Read a fastq file and return the data
    :param fastq_fn:
    :return: seq: (list) a list contains four lines of a fastq file
    """
    try:
        f = open(fastq_fn, "r")
        seq = []
        for line in f:
            seq.append(line)
        f.close()
        return seq
    except IOError:
        print("File not found.")


def getBasename(filename, extention):
    """
    Get basename for a file contains format extension.
    :param filename: (string) filename/path
    :return: basename of a file without path or extension
    """
    filename = filename.split("/")[-1]
    base = filename.replace(extention, "")
    return base


def searchFastq(fast5_fn, fastq_path):
    """
    Retrieve a fastq file base on its correspond fast5 file.
    :param fast5_fn: (string) fast5 file name
    :param fastq_path: (string) path where fastq files are stored
    :return: seq: (list) sequences stores in wanted fastq file
    """
    for fastqfile in os.listdir(fastq_path):
        base = fastqfile.replace(".fastq", "")
        if base == getBasename(fast5_fn, ".fast5"):
            seq = readFastq(fastq_path + fastqfile)
            return seq


def getID(fastq_seq):
    """
    Retrieve sequence ID for a read.
    :param fastq_seq: (list) a list contains four lines of a fastq file
    :return: id: (string) a sequence id
    """
    id = fastq_seq[0].strip().split()[0].replace("@", "")
    return id