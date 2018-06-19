# Author: Haotian Teng
import numpy as np
from src.nanoporereads import *
from src.snps import *
from src.handlefiles import save_objects


### Initialization
AVERAGE_LEN = 4 * 10e3  # Average length of sequence READS
STD_LEN = 10e3  # Standard deviation of the length of sequence READS
SNP_INTERVAL = 3 * 10e2  # Average interval of the SNP among the read
ERROR_RATE = 0.2  # Basecalling error rate
REGION_LEN = 60 * 10e3  # The total length of the sequencing region(Chromsome)
BASE_PROP = [0.25, 0.25, 0.25, 0.25]  # Probability of A,G,C,T
BASE = ['A', 'G', 'C', 'T']  # base list
MATERNAL_SNPS = []  # True maternal SNP list
PATERNAL_SNPS = []  # True paternal SNP list
SNP_POS = []  # True SNP position among the region
READS = []  # The collection of the SNP of READS, it is a list of SNP list,
# e.g. [ [First read SNPS], [Second read SNPS], ...]
LEN_FUNC = np.random.normal
START_FUNC = np.random.uniform
SNP_FUNC = np.random.uniform
SNP_POS_FUNC = np.random.poisson

### True Haplotype
SNP_START = 0
while SNP_START < REGION_LEN:
    add = SNP_POS_FUNC(SNP_INTERVAL)
    SNP_START += add
    SNP_POS.append(SNP_START)
    maternal_roll = SNP_FUNC(low=0, high=np.sum(BASE_PROP))
    maternal_prop = 0
    for idx, prop in enumerate(BASE_PROP):
        maternal_prop += prop
        if maternal_prop >= maternal_roll:
            MATERNAL_SNPS.append(BASE[idx])
            maternal_idx = idx
            break
    paternal_roll = SNP_FUNC(low=0, high=(np.sum(BASE_PROP) - BASE_PROP[maternal_idx]))
    paternal_prop = 0

    for idx, prop in enumerate(BASE_PROP):
        if idx == maternal_idx:
            continue
        else:
            paternal_prop += prop
            if paternal_prop >= paternal_roll:
                PATERNAL_SNPS.append(BASE[idx])
                break
PATERNAL_SNPS = np.asarray(PATERNAL_SNPS)  # 201
MATERNAL_SNPS = np.asarray(MATERNAL_SNPS)  # 201
SNP_POS = np.asarray(SNP_POS)

# Sorting
sort_index = np.argsort(SNP_POS)
PATERNAL_SNPS = PATERNAL_SNPS[sort_index]
MATERNAL_SNPS = MATERNAL_SNPS[sort_index]
SNP_POS = SNP_POS[sort_index]

# Randomly decide (ref-alt)
SNP_OBJ = []
for idx, pos in enumerate(SNP_POS):
    if np.random.rand() >= .5:
        SNP_OBJ.append(SNP(19, idx, pos, PATERNAL_SNPS[idx], MATERNAL_SNPS[idx], "1|0"))
    else:
        SNP_OBJ.append(SNP(19, idx, pos, MATERNAL_SNPS[idx], PATERNAL_SNPS[idx], "0|1"))

# Packing
PM_SNP = np.asarray(list(zip(PATERNAL_SNPS, MATERNAL_SNPS, SNP_POS)),
                    dtype=[('paternal', 'U1'), ('maternal', 'U1'), ('pos', 'i4')])


### Read Simulation
def generate_uniform_transfer_matrix(error_rate):
    transfer_matrix = np.zeros([4, 4])
    for i in range(4):
        transfer_matrix[i] = error_rate / 3
        transfer_matrix[i][i] = 1 - error_rate
    return transfer_matrix


def substitute_base(base_list, transfer_matrix):
    # transfer_matrix is a 4x4 matrix which each element
    # a[_][j] is the probability of _ base transfer to j base.
    bases = ['A', 'G', 'C', 'T']
    transfer_list = []
    for base in base_list:
        base_idx = bases.index(str(base))
        transfer_base = np.random.choice(4, p=transfer_matrix[base_idx])
        transfer_list.append(bases[transfer_base])
    return transfer_list


READ_START = 0
READ_LEN = 0
READS_NUM = 1000
READ_REGION = [2 * 10e3, 10 * 10e3]  # The sub region of the READS, the read is restricted in this region
error_matrix = generate_uniform_transfer_matrix(ERROR_RATE)

ID = 0
R_OBJ = []
READS = []

while True:
    READ_START = int(START_FUNC(low=READ_REGION[0], high=READ_REGION[1]))
    READ_LEN = int(LEN_FUNC(loc=AVERAGE_LEN, scale=STD_LEN))  # 20000 to 50000
    if (READ_START + READ_LEN) > READ_REGION[1]:
        continue
    snp_index = np.where(np.logical_and(PM_SNP['pos'] > READ_START, PM_SNP['pos'] < (READ_START + READ_LEN)))

    p_m = np.random.randint(low=0, high=2)
    if p_m == 0:
        snp_list = PM_SNP[snp_index]['paternal']
    else:
        snp_list = PM_SNP[snp_index]['maternal']
    t_list = substitute_base(snp_list, error_matrix)
    packing = (t_list, snp_index)  # snp_idx: tuple
    READS.append(packing)

    # snp_id: read base. dict
    BASES = {}
    for index, b in enumerate(t_list):
        BASES[snp_index[0][index]] = b

    print(BASES)
    # read objects
    r_obj = NanoporeRead(ID, 19, READ_START, (READ_LEN+READ_START), 60)
    r_obj.set_bases(BASES)
    r_obj.detect_snps(SNP_OBJ)
    R_OBJ.append(r_obj)
    ID += 1

    if len(READS) >= READS_NUM:
        break

save_objects("data/dr1.obj", R_OBJ)
save_objects("data/ds1.obj", SNP_OBJ)

#s = []
#for snp in SNP_OBJ:
    #for read in R_OBJ:
        #if read.start <= snp.pos <= read.end:
            #if snp not in s:
                #s.append(snp)
#print(len(s))  # 27
