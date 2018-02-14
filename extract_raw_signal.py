"""
__author__ = Yao LI
__email__ = liyaoo1012@163.com
__date__ = 07/02/2018
14/02/18
"""
from h5utils import *
import os

f = open("raw_signal.txt", "w")

for filename in os.listdir("/shares/coin/yao.li/data/basecall_pass/"):
    if filename.endswith(".fast5"):
        seg_raw, seg_fastq = get_raw_segment(filename, 100, 200)
        f.write(seg_raw)
        f.write(seg_fastq)
f.close()