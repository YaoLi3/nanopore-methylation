#! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 20/02/2018
"""
import numpy as np


#########################
# Dynamic Time Wrapping #
#########################
def raw_signal_clustering(raw_signals):
    """
    Cluster nanopore sequencing raw signal into two clusters:
    methylation or non-methylation groups
    """
    CLUSTERS = ["Methylation", "Non-methylation"]
    self.data = raw_signals