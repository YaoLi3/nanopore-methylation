import ggplot as gp

from matplotlib import pyplot as plt
from matplotlib import style

from dna_features_viewer import GraphicFeature, GraphicRecord


def models_llhd(pm_llhd):
    """
    Tracking the total likelihood of reads in a model(cluster).
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
    s1, s2 = model.read_results()
    record = GraphicRecord(sequence=ref_seq, sequence_length=len(ref_seq),
                           features=[GraphicFeature(start=0, end=len(s1), strand=+1, color='#ffcccc'),
                                     GraphicFeature(start=0, end=len(s2), strand=+1, color='#ffcccc')])  # change color
    ax, _ = record.plot(figure_width=5)
    record.plot_sequence(ax)
    record.plot_translation(ax, (8, 23), fontdict={'weight': 'bold'})
    ax.figure.savefig('haplotypes.png', bbox_inches='tight')
