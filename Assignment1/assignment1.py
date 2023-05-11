#!/usr/bin/env python3

"""
python3 assignment1.py -n <aantal_cpus> [OPTIONEEL: -o <output csv file>] fastabestand1.fastq [fastabestand2.fastq ... fastabestandN.fastq]
"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.1"

import re
# IMPORTS
import sys
import argparse as ap
import re

# CLASSES
class QualityController:
    """ Desc """
    def __init__(self):
        self.fastqfile = '/commons/Themas/Thema12/HPC/rnaseq_selection.fastq'

    def arg_parse(self):
        argparser = ap.ArgumentParser(description="Script voor Opdracht 1 van Big Data Computing")
        argparser.add_argument("-n", action="store",
                               dest="n", required=True, type=int,
                               help="Aantal cores om te gebruiken.")
        argparser.add_argument("-o", action="store", dest="csvfile", type=ap.FileType('w', encoding='UTF-8'),
                               required=False,
                               help="CSV file om de output in op te slaan. Default is output naar terminal STDOUT")
        argparser.add_argument("fastq_files", action="store", type=ap.FileType('r'), nargs='+',
                               help="Minstens 1 Illumina Fastq Format file om te verwerken")
        return argparser.parse_args()

    def file_handler(self):
        """ Iterate through fastqfile. NB: chunk is assumed to be multiple of 4!"""
        with open(self.fastqfile, 'r') as fastq:
            for line in fastq:
                if re.match(r'^[ACTGN]+$', line):
                    print(line + 'nucleo')
                elif re.match(r'\+', line):
                    print(line + '+')
                else:
                    print("header")

        # with open(fastqfile, 'r') as fastq:
        #     # fastforward
        #     i = 0
        #     while i < start:
        #         fastq.readline()
        #         i += 1
        #
        #     results = []
        #     counter = 0
        #     while counter < chunk:
        #         header = fastq.readline()
        #         nucleotides = fastq.readline()
        #         strand = fastq.readline()
        #         qual = fastq.readline()
        #         counter += 4
        #
        #         if not (qual):
        #             # we reached the end of the file
        #             break
        #         for j, c in enumerate(qual):
        #
        #             try:
        #                 results[j] += ord(c) - 33
        #             except IndexError:
        #                 results.append(ord(c) - 33)
        #
        #     return [(phredscore / (counter / 4)) for phredscore in results]


# FUNCTIONS
def y():
    """ Desc """
    return 0


# MAIN
def main():
    """ Main function """
    # FINISH
    qc = QualityController()
    fastq = "/commons/Themas/Thema12/HPC/rnaseq_selection.fastq"
    qc.file_handler()
    return 0


if __name__ == '__main__':
    main()
