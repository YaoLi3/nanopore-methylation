import ggplot as gp
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import style
from dna_features_viewer import GraphicFeature, GraphicRecord
from src.handlefiles import read_fastq

CHR19_REF_SEQ = read_fastq("../data/GRCh38_full_analysis_set_plus_decoy_hla.fa")


def models_llhd(pm_llhd):
    """
    Tracking the total likelihood of READS in a model(cluster).
    :param pm_llhd: (np.array) matrix stores read likelihood in every model/cluster.
    :param type (

    x axis: iteration time
    y axis: sum likelihood log value
    """
    p = gp.ggplot(gp.aes(x="iteration num", y="log value"), data=pm_llhd)\
        +gp.geom_point(color="blue")\
        +gp.ggtitle(u"model likelihood")
    print(p)


def save_scatter_plot_fig(content, fn):
    """
    :param content: (np.array)
    :param fn: (str) file name
    """
    style.use("ggplot")
    plt.plot(content)
    plt.savefig(fn, bbox_inches='tight')
    plt.show()


def haplotype_blocks_fig(model, ref_seq):
    s1, s2 = model.align_alleles()
    record = GraphicRecord(sequence=ref_seq, sequence_length=len(ref_seq),
                           features=[GraphicFeature(start=0, end=len(s1), strand=+1, color='#ffcccc'),
                                     GraphicFeature(start=0, end=len(s2), strand=+1, color='#cffccc')])
    ax, _ = record.plot(figure_width=5)
    record.plot_sequence(ax)
    record.plot_translation(ax, (8, 23), fontdict={'weight': 'bold'})
    ax.figure.savefig('haplotypes.png', bbox_inches='tight')


def plot_sr_dict(sr_dict, fn):
    """
    x: snp_id
    y: mapped reads numbers
    """
    snps = sorted(list(sr_dict.keys()))
    l1 = []
    l2 = []
    for snp in sr_dict:
        l1.append(sr_dict[snp][0])
        l2.append(sr_dict[snp][1])
    #my_x_ticks = np.arange(9490, 50000, 20)
    #plt.xticks(my_x_ticks)
    plt.scatter(snps, l1, label="1st haplotype")
    plt.scatter(snps, l2, c="red", label="2nd haplotype")
    plt.xlabel(u"SNP id")
    plt.ylabel(u"Reads num")
    plt.legend(loc='upper right')
    plt.savefig(fn, bbox_inches="tight")
    plt.show()
