#! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 20/02/2018
"""
import os
import h5py
import numpy as np


def extract_fastq(fast5_fn):
    """
    Extract fastq sequence from a given fast5 file,
    and save the fastq sequence in a file with the same name.
    :param fast5_fn: (string) fast5 file name
    """
    try:
        f = h5py.File(fast5_fn, "r")
        seq = f['Analyses']['Basecall_1D_001']['BaseCalled_template']['Fastq'].value
        if seq != "":
            fastq_file = open("%s.fastq"% fast5_fn.replace(".fast5", ""), "wb")
            fastq_file.write(seq)
            fastq_file.close()
    except IOError:
        print("File not found.")


def get_fast5_id(fast5_fn):
    """
    :param fast5_fn:
    :return:
    """
    try:
        root = h5py.File(fast5_fn, "r")
        atr = root["/Analyses/Basecall_1D_000/BaseCalled_template/"].value
        return atr
    except IOError:
        print("File not found.")


def read_fastq(fastq_fn):
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


def get_basename(filename, extention):
    """
    Get basename for a file contains format extension.
    :param filename: (string) filename/path
    :return: basename of a file without path or extension
    """
    if extention in filename:
        filename = filename.split("/")[-1]
        base = filename.replace(extention, "")
        return base
    else:
        raise ValueError("Wrong file type.")


def search_fastq(fast5_fn, fastq_path):
    """
    Retrieve a fastq file base on its correspond fast5 file.
    :param fast5_fn: (string) fast5 file name
    :param fastq_path: (string) path where fastq files are stored
    :return: seq: (list) sequences stores in wanted fastq file
    """
    try:
        for fastqfile in os.listdir(fastq_path):
            base = fastqfile.replace(".fastq", "")
            if base == get_basename(fast5_fn, ".fast5"):
                seq = read_fastq(fastq_path + fastqfile)
                return seq
    except IOError:
        raise IOError("{} not exist.".format(fastq_path))


def get_id(fastq_seq):
    """
    Retrieve sequence ID for a read.
    :param fastq_seq: (list) a list contains four lines of a fastq file
    :return: id: (string) a sequence id
    """
    if "@" in fastq_seq[0]:
        id = fastq_seq[0].strip().split()[0].replace("@", "")
        return id
    else:
        raise ValueError("Wrong fastq format.")


def getID(fast5_f):
    """
    @author: haotian.teng
    :param fast5_f:
    """
    with h5py.File(fast5_f, 'r') as handle:
        protocal_id = handle['UniqueGlobalKey/tracking_id'].attrs['protocol_run_id']
    return protocal_id


def load_npy(npy_fn):
    """
    :param npy_fn:
    :return:
    """
    return np.load(npy_fn)[0]