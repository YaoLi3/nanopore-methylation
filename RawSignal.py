"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 07/02/2018
"""
import os
from h5utils import *
from imprinted_reads_less import *
from handleFiles import *

if __name__ == "__main__":
    DATA = NanoporeReads("/shares/coin/yao.li/minimap2/merged.sam", "19")
    DATA.getReads()  # 45946 reads
    o = DATA.findImprinted(ImprintedRegions
                           ("/shares/coin/yao.li/data/mart_export.txt").getRegions(),
                           0, True, "find_imprinted_result.txt")


    f = open("results.txt", "w")
    # Extract raw signal from fast5 files
    for fst5 in os.listdir("/shares/coin/yao.li/data/basecall_pass/"):
        if fst5.endswith(".fast5"):
            try:
                qpath = "/shares/coin/yao.li/data/fastq/"
                poses = o[getID(searchFastq(fst5, qpath))]
                raw, fastq = get_raw_segment(fst5, poses[0], poses[1])
                np.save("/shares/coin/yao.li/raw_signal/{}.npy".format(fst5.replace(".fast5", "")), (raw, fastq))
            except KeyError:
                f.write("this file does not have basecalled_template")
                continue
            except TypeError:
                f.write("Read not found in imprinted regions.")
                continue
            except StopIteration:
                f.write("This is the end.")
    f.close()