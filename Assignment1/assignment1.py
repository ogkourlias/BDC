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
import sys
import multiprocessing as mp
import pandas as pd
import numpy as np


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
    Class containg functiosn to calculate the average Phred scores for positions.
    """

    def __init__(self, args, chunk_size=50000):
        self.files = args.fastq_files
        self.chunk_size = chunk_size
        self.cores = args.n
        self.fastq_files = args.fastq_files
        self.csvfile = args.csvfile

    def calculate(self):
        """
        Function containing the primary calculation pipeline for  given files.
        """
        for file in self.files:
            with file as fastq:
                self.line_walker(fastq)

    def line_walker(self, file):
        """
        Walks through the lines the given input
        :param file:
        """
        with open(file) as fastq:
            pool = mp.Pool(self.cores)
            score_lines = []
            lines = fastq.readlines()
            for i, line in enumerate(lines):
                if (i + 1) % 4 == 0:
                    score_lines.append(line.strip())

            res = pool.map(self.score_getter, [score_lines[i:i + self.chunk_size]
                                               for i in range(0, len(score_lines), self.chunk_size)])

            score_list = [np.array(result[0]) for result in res]
            pos_list = [np.array(result[1]) for result in res]
            arr1 = np.array(
                [np.pad(row, (0, len(max(score_list, key=len)) - len(row)), 'constant') for row in score_list])
            arr2 = np.array(
                [np.pad(row, (0, len(max(score_list, key=len)) - len(row)), 'constant') for row in pos_list])
            scores_summed = np.sum(arr1, axis=0).tolist()
            pos_summed = np.sum(arr2, axis=0).tolist()
            output = [score_sum / pos_sum for score_sum, pos_sum in zip(scores_summed, pos_summed)]
            self.write_output(output)

    def score_getter(self, lines):
        """
        Gets the scores and indexes for lines
        :param lines:
        :return:
        """
        pos_scores = []
        pos_counts = []
        pos_dict = {}
        for line in lines:
            for i in range(len(line)):
                if len(pos_scores) > i:
                    pos_dict[i] = pos_dict[i] + 1
                    pos_counts[i] += 1
                    pos_scores[i] += (ord(line[i]) - 33)
                else:
                    pos_dict[i] = 1
                    pos_counts.append(1)
                    pos_scores.append(ord(line[i]) - 33)

        return pos_scores, pos_counts

    def write_output(self, output):
        """Write scores to csv file."""
        df = pd.DataFrame(output)
        if self.csvfile:
            df.to_csv(self.csvfile, header=False)
        else:
            df.to_csv(sys.stdout, header=False)


# MAIN
def main():
    """Main function"""
    args = arg_parse()
    calc = AvgCalc(args)
    calc.calculate()


if __name__ == "__main__":
    main()
