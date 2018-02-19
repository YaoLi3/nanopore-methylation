"""
__author__ = Yao LI
__email__ = liyaoo1012@163.com
__date__ = 07/02/2018
14/02/18
"""
import os
from h5utils import *
from imprinted_reads import *


def imprinted_raw(fast5_file):
    # get the read id
    rid = extract_fastq(fast5_file)[0].strip().split()[0]
    try:
        # the read is imprinted
        gene_name, status, read_pos, ref_pos = o[rid]
        # extract raw sigal
        seg_raw, seg_fastq = get_raw_segment(fast5_file, 100, 200)
        raw = seg_raw[read_pos[0]: read_pos[1]]
        return raw, seg_fastq
    except KeyError:  # if the read does not overlapped with any human imprinted region
        print("key error, read does not qualified")
        return


if __name__ == "__main__":
    # Get human imprinted regions
    IR = ImprintedRegions("/home/yaoli/Desktop/ppth/data/mart_export.txt")
    imprinted_regions = IR.getRegions()

    # Retrieve Nanopore reads
    DATA = NanoporeReads("/home/yaoli/Desktop/ppth/data/merged.sam")
    DATA.getReads()  # 45946 reads
    o = DATA.findImprinted(imprinted_regions, 0, True, "find_imprinted_result.txt")


    # Extract raw signal from fast5 files
    for filename in os.listdir("/shares/coin/yao.li/data/basecall_pass/"):
        if filename.endswith(".fast5"):
            try:
                name = "/shares/coin/yao.li/data/basecall_pass/" + filename
                raw, fastq = imprinted_raw(name)
                np.save("/shares/coin/yao.li/raw_signal/" + filename, (raw, fastq))
            except KeyError:
                print("this file does not have basecalled_template")
                continue
            except StopIteration:
                print("THis is the end.")


