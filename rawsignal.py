#! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 07/02/2018
"""
from h5utils import *
from nanoporereads import *
from handlefiles import *


def get_raw_dirc(directory, savepath, ir_pos):
    """
    :param directory: (string) the folder of fast5 files
    :param qpath: (string) path of fastq files
    :param savepath: (string) path of output numpy files
    :param ir_pos: (dictionary) a NanoporeReads overlap dict
    :return: (dict) raw signal and its fastq sequence
    """
    raw_signal = {}
    # Extract raw signal from fast5 files
    for fst5 in os.listdir(directory):
        try:
            if fst5.endswith(".fast5"):
                sid = getID(directory+fst5)
                if sid in ir_pos:
                    poses = ir_pos[sid][2]
                    raw, fastq = get_raw_segment(directory+fst5, poses[0], poses[1], basecall_group='Basecall_1D_001')
                    raw_signal[fst5] = (raw, fastq)
                    np.save(savepath + fst5.replace(".fast5", ".npy"), (raw, fastq))
        except StopIteration:
            print("This is the end.")
    return raw_signal


def find_haplotype(raw_signals):
    """For each raw_signal array, decide its haplotype."""
    #TODO: use read id should be straightforward. need a class?
    pass
