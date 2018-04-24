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
            pickle.dump(read, f, protocol=2)


def load_pickle_file(filename):
    """Read a pickle file, return a generator"""
    with open(filename, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break


def load_objects(filename):
    """Return a list of objects"""
    return list(load_pickle_file(filename))


def compare_lists(list1, list2):
    """Compare to list1, changes happened in list2.
    Measure the difference between 2 lists."""
    add = []
    miss = []
    for ele in list1:
        if ele not in list2:
            miss.append(ele)
    for ele in list2:
        if ele not in list1:
            add.append(ele)

    por = len(miss)/len(list1)
    print("{}% of 1st list elements are missing in the 2nd list.".format(por))
    print("{} elements in 2nd list are additions to the 1st list.".format(len(add)))
    return add, miss


"""
Created on Tue Feb 13 16:04:52 2018

@author: haotianteng
"""


def get_raw_segment(fast5_fn, start_base_idx, end_base_idx, basecall_group='Basecall_1D_000',
                    basecall_subgroup='BaseCalled_template'):
    """
    Get the raw signal segment given the start and end snp_id of the sequence.
    fast5_fn: input fast5 file name.
    start_base_idx: start snp_id of the sequence (0-based)
    end_base_idx: end snp_id of the sequence (the snp_id is included)
    basecall_group: group name to search for base information.
    basecall_subgroup: sub grou#!p name to search for base information.

    e.g.
        get_raw_segment('test.fast5', 0, 10)
        Will return the signal corresponded to the 0-10 bases(The 0th and 10th base are both included.)


    """
    with h5py.File(fast5_fn, 'r') as root:
        base = root['Analyses/{}/BaseCalled_template'.format(basecall_group)]
        fastq = base['Fastq'].value.split()[2]
        seg = fastq[start_base_idx:end_base_idx]
        event_h = base['Events']
        events = event_h.value
        raw_h = list(root['/Raw/Reads'].values())
        raw = raw_h[0]['Signal']
        start_time = None
        if (type(events[0][1]) is np.float64) or (type(events[0][1]) is np.float32):
            start_time = event_h.attrs['start_time']
        pos = list()
        pos_idx = 0
        for event in events:
            pos_idx += event[5]
            pos.append(pos_idx)
        start_idx = next(x[0] for x in enumerate(pos) if x[1] >= start_base_idx)
        end_idx = next(x[0] - 1 for x in enumerate(pos) if x[1] > end_base_idx)
        if start_time is None:
            raw_start = events[start_idx][1]
            raw_end = events[end_idx][1]
        else:
            raw_start = int((events[start_idx][1] - start_time) / 0.00025)
            raw_end = int((events[end_idx][1] - start_time) / 0.00025)
        seg_raw = raw[raw_start:raw_end]
    return seg_raw, seg