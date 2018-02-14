#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 16:04:52 2018

@author: haotianteng
"""

import h5py
import numpy as np

def get_raw_segment(fast5_fn,start_base_idx, end_base_idx, basecall_group = 'Basecall_1D_000',basecall_subgroup = 'BaseCalled_template'):
    """
    Get the raw signal segment given the start and end index of the sequence.
    fast5_fn: input fast5 file name.
    start_base_idx: start index of the sequence (0-based)
    end_base_idx: end index of the sequence (the index is included)
    basecall_group: group name to search for base information.
    basecall_subgroup: sub group name to search for base information.
    
    e.g.
        get_raw_segment('test.fast5', 0, 10)
        Will return the signal corresponded to the 0-10 bases(The 0th and 10th base are both included.)
        
    
    """
    with h5py.File(fast5_fn,'r') as root:
        base = root['Analyses/Basecall_1D_000/BaseCalled_template']
        fastq = base['Fastq'].value.split()[2]
        seg =fastq[start_base_idx:end_base_idx]    
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
        end_idx = next(x[0]-1 for x in enumerate(pos) if x[1] > end_base_idx)
        if start_time is None:
            raw_start = events[start_idx][1]
            raw_end = events[end_idx][1]
        else:
            raw_start = int((events[start_idx][1]-start_time)/0.00025)
            raw_end = int((events[end_idx][1]-start_time)/0.00025)
        seg_raw = raw[raw_start:raw_end]
    return seg_raw, seg
        
if __name__=="__main__":
    
    seg_raw,seg_fastq = get_raw_segment('/home/yaoli/PycharmProjects/project/data/PLSP61583_20160920_FNFAB390088_MN17048_sequencing_run_Hum_94_62579_ch268_read765_strand.fast5',100,200)
    print(seg_raw)
    print(seg_fastq)