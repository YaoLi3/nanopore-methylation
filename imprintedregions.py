#! /usr/bin/env python
"""
__author__ = Yao LI
__email__ = yao.li.binf@gmail.com
__date__ = 08/02/2018
"""


#############################
#  Human Imprinted Regions  #
#############################
class ImprintedRegions:
    """
    Define the imprinted regions on human genome (EnSembl hg38).
    Save the data in a dictionary.
    """

    def __init__(self, filename):
        """
        :param filename: (string) path of the file stores the
                         information of certain imprinted regions.
        """
        self.regions = {}
        file = open(filename, "r")
        for line in file:
            line = line.strip().split("\t")
            self.regions[line[0]] = line[1:]
        file.close()

    def get_regions(self):
        """
        :return: (dictionary) key = gene names
        """
        return self.regions