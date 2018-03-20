import os
import h5py


def extract_fastq(name):
    f = h5py.File("/shares/coin/yao.li/test_data/%s.fast5" % name, "r")
    seq = f['Analyses']['Basecall_1D_001']['BaseCalled_template']['Fastq'].value
    if seq != "":
        fastq_file = open("%s.fastq" % name, "wb")
        fastq_file.write(seq)
        fastq_file.close()


if __name__ == "__main__":
    dirt = os.listdir("/shares/coin/yao.li/test_data/")
    for reads in dirt:
        reads_name = reads.replace(".fast5", "")
        extract_fastq(reads_name)