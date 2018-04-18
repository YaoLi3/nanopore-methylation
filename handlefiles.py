#! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 20/02/2018
"""
import os
import h5py
import numpy as np
import pickle


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
    :return: atr
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
    Retrieve sequence ID for a READs.
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
    :return: (numpy) raw signal array
    """
    return np.load(npy_fn)[0]


def read_find_ip_results(fn):
    """
    Read a text file where save ... data, tab delimiter ...
    :param fn: result file name
    """
    f = open(fn, "r")
    data = {}
    for line in f:
        line = line.strip().split("\t")
        if len(line) == 7:
            id, gene, info, read_pos, ref_pos, thrhld, seq = line
            read_pos = read_pos.strip("()").split(",")
            ref_pos = ref_pos.strip("()").split(",")
            data[id] = (gene, info, read_pos, ref_pos, thrhld, seq)


def get_positions(file_n):
    """
    :param file_n: (string) file name
    :return: (dict) key = reads id, value = (start, end) genome positions of the READs
    """
    poses = {}
    f = open(file_n, "r")
    for line in f:
        line = line.strip().split("\t")
        if len(line) == 7:
            pos = line[4].strip("()").split(",")
            poses[line[0]] = (int(pos[0]), int(pos[1][1:]))
    return poses


def find_snps_in_read(read_f, snp_data):
    """
    :param read_f: (string) find_imprinted_result file name
    :param snp_data: (list) [(SNPs genome position, a status, b status)]
    :return: (dict) data: key = id, value = info. key = SNPs, seq.
    """
    data = {}
    ip_reads_regions = get_positions(read_f)
    for id in ip_reads_regions:
        snps = []
        start, end = ip_reads_regions[id]
        for snp in snp_data:
            snp_p = snp_data[snp][1]
            if start <= snp_p <= end:
                snps.append(snp_data[snp])
        if not snps == []:
            data[id] = snps
    return data


def save_objects(filename, reads):
    """Save NanoporeRead objects into a file"""
    with open(filename, "wb") as f:
        for read in reads:
            pickle.dump(read, f)


def load_objects_file(filename):
    """Read a obj file, return a list of NanoporeRead objects"""
    with open(filename, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break
