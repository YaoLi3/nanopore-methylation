"""
__author__ = Yao LI
__email__ = liyaoo1012@163.com
__date__ = 07/02/2018
14/02/18
"""
import sys
import os
from h5utils import *
from imprinted_reads_less import *


def imprinted_raw(fast5):
    # get the read id
    rid = extract_fastq(fast5)[0].strip().split()[0]
    try:
        # the read is imprinted
        gene_name, status, read_pos, ref_pos = o[rid]
        # extract raw signal
        seg_raw, seg_fastq = get_raw_segment(fast5, read_pos[0], read_pos[1])
        return seg_raw, seg_fastq
    except KeyError:  # if the read does not overlapped with any human imprinted region
        print("key error, read does not qualified")
    except IndexError:
        print("Raw signal does not have adequate amount of data.")


if __name__ == "__main__":
    # Get human imprinted regions
    IR = ImprintedRegions("/shares/coin/yao.li/data/mart_export.txt")
    imprinted_regions = IR.getRegions()

    # Retrieve Nanopore reads
    DATA = NanoporeReads("/shares/coin/yao.li/minimap2/merged.sam", "19")
    DATA.getReads()  # 45946 reads
    o = DATA.findImprinted(imprinted_regions, 0, True, "find_imprinted_result.txt")

    # Extract raw signal from fast5 files
    for basename in os.listdir("/shares/coin/yao.li/data/basecall_pass/"):
        if basename.endswith(".fast5"):
            try:
                path = "/shares/coin/yao.li/data/basecall_pass/"
                raw, fastq = imprinted_raw(path+basename)
                np.save("/shares/coin/yao.li/raw_signal/{}.npy".format(basename.replace(".fast5")), (raw, fastq))
            except KeyError:
                print("this file does not have basecalled_template")
                continue
            except TypeError:
                print("Read not found in imprinted regions.")
                continue
            except StopIteration:
                print("This is the end.")
