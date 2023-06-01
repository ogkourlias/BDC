#!/usr/bin/env python3

"""python3 assignment1.py -n <aantal_cpus> [OPTIONEEL: -o <output csv file>] fastabestand1.fastq
[fastabestand2.fastq ... fastabestandN.fastq]"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.1"

# IMPORTS
import argparse as ap
from pathlib import Path
import csv
import multiprocessing as mp


# FUNCTIONS
def arg_parse():
    """Parse arguments"""
    argparser = ap.ArgumentParser(description="Script voor Opdracht 1 van Big Data Computing")
    argparser.add_argument("-n", action="store",
                           dest="n", required=True, type=int,
                           help="Aantal cores om te gebruiken.")
    argparser.add_argument("-o", action="store", dest="csvfile", type=Path,
                           required=False,
                           help="CSV file om de output in op te slaan. Default is output naar terminal STDOUT")
    argparser.add_argument("fastq_files", action="store", type=Path, nargs='+',
                           help="Minstens 1 Illumina Fastq Format file om te verwerken")
    return argparser.parse_args()


class AvgCalc:
    """
    t
    """

    def __init__(self, files, chunk_size=2, cores=12):
        self.files = files
        self.chunk_size = chunk_size
        self.cores = cores

    def receiver(self):
        """
        :return:
        """
        return 0

    def calculate(self):
        for file in self.files:
            with file as fastq:
                self.line_walker(fastq)

    def line_walker(self, file):
        with open(file) as fastq:
            pool = mp.Pool(self.cores)
            score_lines = []
            lines = fastq.readlines()
            for i, line in enumerate(lines):
                if (i + 1) % 4 == 0:
                    score_lines.append(line.strip())

            scores = pool.map(self.score_getter, [score_lines[i:i + self.chunk_size]
                                                  for i in range(0, len(score_lines), self.chunk_size)])

            print(scores)


    def score_getter(self, lines):
        scores = []
        for line in lines:
            scores = [ord(char) - 33 for char in line]
        return scores


# MAIN
def main():
    """Main function"""
    args = arg_parse()
    calc = AvgCalc(args.fastq_files)
    calc.calculate()


if __name__ == "__main__":
    main()
