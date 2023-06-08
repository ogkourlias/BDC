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
from demo import *
import sys
import multiprocessing as mp
import pandas as pd
import numpy as np


# FUNCTIONS
def arg_parse():
    argparser = ap.ArgumentParser(
        description="Script voor Opdracht 2 van Big Data Computing;  Calculate PHRED scores over the network.")
    mode = argparser.add_mutually_exclusive_group(required=True)
    mode.add_argument("-s", action="store_true", help="Run the program in Server mode; see extra options needed below")
    mode.add_argument("-c", action="store_true", help="Run the program in Client mode; see extra options needed below")
    server_args = argparser.add_argument_group(title="Arguments when run in server mode")
    server_args.add_argument("-o", action="store", dest="csvfile", type=Path,
                             required=False,
                             help="CSV file om de output in op te slaan. Default is output naar terminal STDOUT")
    server_args.add_argument("fastq_files", action="store", type=Path, nargs='*',
                             help="Minstens 1 Illumina Fastq Format file om te verwerken")
    server_args.add_argument("--chunks", action="store", type=int, required=True)

    client_args = argparser.add_argument_group(title="Arguments when run in client mode")
    client_args.add_argument("-n", action="store",
                             dest="n", required=False, type=int,
                             help="Aantal cores om te gebruiken per host.")
    client_args.add_argument("--host", action="store", type=str, help="The hostname where the Server is listening")
    client_args.add_argument("--port", action="store", type=int, help="The port on which the Server is listening")
    return argparser.parse_args()

class AvgCalc:
    """
    Class containg functiosn to calculate the average Phred scores for positions.
    """
    def calculate(self, files):
        """
        Function containing the primary calculation pipeline for  given files.
        """
        for file in files:
            with file as fastq:
                self.line_walker(fastq)

    def line_walker(self, file, cores = 12, chunk = 5000):
        """
        Walks through the lines the given input
        :param file:
        """
        with open(file) as fastq:
            pool = mp.Pool(cores)
            score_lines = []
            lines = fastq.readlines()
            for i, line in enumerate(lines):
                if (i + 1) % 4 == 0:
                    score_lines.append(line.strip())

            res = pool.map(self.score_getter, [score_lines[i:i + chunk]
                                               for i in range(0, len(score_lines), chunk)])

            score_list = [np.array(result[0]) for result in res]
            pos_list = [np.array(result[1]) for result in res]
            arr1 = np.array(
                [np.pad(row, (0, len(max(score_list, key=len)) - len(row)), 'constant') for row in score_list])
            arr2 = np.array(
                [np.pad(row, (0, len(max(score_list, key=len)) - len(row)), 'constant') for row in pos_list])
            scores_summed = np.sum(arr1, axis=0).tolist()
            pos_summed = np.sum(arr2, axis=0).tolist()
            output = [score_sum / pos_sum for score_sum, pos_sum in zip(scores_summed, pos_summed)]
            self.write_output(output, "test.csv")

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

    def calc_final(self, args, res):
        score_list = [np.array(result[0]) for result in res]
        pos_list = [np.array(result[1]) for result in res]
        arr1 = np.array(
            [np.pad(row, (0, len(max(score_list, key=len)) - len(row)), 'constant') for row in score_list])
        arr2 = np.array(
            [np.pad(row, (0, len(max(score_list, key=len)) - len(row)), 'constant') for row in pos_list])
        scores_summed = np.sum(arr1, axis=0).tolist()
        pos_summed = np.sum(arr2, axis=0).tolist()
        output = [score_sum / pos_sum for score_sum, pos_sum in zip(scores_summed, pos_summed)]
        self.write_output(output, "test1.csv")

    def write_output(self, output, csvfile):
        """Write scores to csv file."""
        df = pd.DataFrame(output)
        if csvfile:
            df.to_csv(csvfile, header=False)
        else:
            df.to_csv(sys.stdout, header=False)

# MAIN
def main():
    """Main function"""
    args = arg_parse()
    avgCalc = AvgCalc()

    server = mp.Process(target=runserver, args=(avgCalc.line_walker, args.fastq_files))
    server.start()
    time.sleep(1)
    client = mp.Process(target=runclient, args=(4,))
    client.start()
    server.join()
    client.join()


if __name__ == "__main__":
    main()
