#!/usr/bin/env python3

"""
Script to calculate PHRED scores over the network in either a client or server mode. Outputs a CSV file with the scores.
"""

__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.3"

import argparse as ap
import csv
from pathlib import Path
import sys
import pandas as pd
import numpy as np


def arg_parse():
    """
    Parse the command line arguments.
    :return: Parsed command line arguments.
    """
    argparser = ap.ArgumentParser(
        description="Script to calculate PHRED scores over the network.")

    mode = argparser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--chunker", action="store_true",
                      help="Run the program in Server mode; see extra options needed below")
    mode.add_argument("--combine", action="store_true",
                      help="Run the program in Client mode; see extra options needed below")

    server_args = argparser.add_argument_group(title="Arguments when run in server mode")
    server_args.add_argument("-o", action="store", dest="csvfile", type=Path, required=False,
                             help="CSV file to save the output. Default is output to terminal STDOUT")
    server_args.add_argument("fastq_files", action="store", type=str, nargs='*',
                             help="At least 1 Illumina Fastq Format file to process")
    server_args.add_argument("--chunks", action="store", type=int, required=False,
                             help="Chunk size for processing")

    client_args = argparser.add_argument_group(title="Arguments when run in client mode")
    client_args.add_argument("-n", action="store", dest="n", required=False, type=int,
                             help="Number of cores to use per host.")
    client_args.add_argument("--host", action="store", type=str,
                             help="The hostname where the Server is listening")
    client_args.add_argument("--port", action="store", type=int,
                             help="The port on which the Server is listening")

    return argparser.parse_args()


class AvgCalc:
    """
    Class containing functions to calculate the average Phred scores for positions.
    """

    def __init__(self, args):
        self.chunk_size = 1000
        self.files = None
        self.csvfile = None

    def work_division(self, lines):
        """
        Divide the lines into chunks of self.chunk_size.
        :param lines: Lines to be divided.
        :return: List of chunks.
        """
        return [lines[i:i + self.chunk_size] for i in range(0, len(lines), self.chunk_size)]

    def walk_through_files(self, files=None):
        """
        Process the provided files.
        :param files: List of file paths. If None, uses self.files.
        :return: List of outputs for each file.
        """
        files = files or self.files
        return [self.walk_through_lines(line_list) for line_list in self.files_handler(files)]

    def walk_through_lines(self, line_list):
        """
        Process the given lines and return the average scores.
        :param line_list: Lines to process.
        :return: Average scores.
        """
        work = self.work_division(line_list)
        res = [self.score_getter(chunk) for chunk in work]
        score_list = [np.array(result[0]) for result in res]
        pos_list = [np.array(result[1]) for result in res]
        max_len = len(max(score_list, key=len))

        arr1 = np.array([np.pad(row, (0, max_len - len(row)), 'constant') for row in score_list])
        arr2 = np.array([np.pad(row, (0, max_len - len(row)), 'constant') for row in pos_list])

        scores_summed = np.sum(arr1, axis=0).tolist()
        pos_summed = np.sum(arr2, axis=0).tolist()

        output = [score_sum / pos_sum for score_sum, pos_sum in zip(scores_summed, pos_summed)]
        return output

    def score_getter(self, lines):
        """
        Get the scores and indexes for lines.
        :param lines: Lines to get scores for.
        :return: Tuple of scores and counts.
        """
        pos_scores = []
        pos_counts = []
        for line in lines:
            for i in range(len(line)):
                if i < len(pos_scores):
                    pos_counts[i] += 1
                    pos_scores[i] += (ord(line[i]) - 33)
                else:
                    pos_counts.append(1)
                    pos_scores.append(ord(line[i]) - 33)

        return pos_scores, pos_counts

    def files_handler(self, files):
        """
        Process the given files and return lines for each.
        :param files: List of file paths.
        :return: List of lines for each file.
        """
        return [[line.strip() for i, line in enumerate(open(file).readlines()) if (i + 1) % 4 == 0] for file in files]

    def write_output(self, output):
        """
        Write scores to a csv file.
        :param output: Output to write.
        """
        df = pd.DataFrame(output)
        df.to_csv(self.csvfile or sys.stdout, header=False)

    def chunk_handler(self, chunk):
        return self.score_getter(chunk.readlines())

    def chunk_combine(self, res):
        score_list = [np.array(result[0]) for result in res]
        pos_list = [np.array(result[1]) for result in res]
        max_len = len(max(score_list, key=len))

        arr1 = np.array([np.pad(row, (0, max_len - len(row)), 'constant') for row in score_list])
        arr2 = np.array([np.pad(row, (0, max_len - len(row)), 'constant') for row in pos_list])

        scores_summed = np.sum(arr1, axis=0).tolist()
        pos_summed = np.sum(arr2, axis=0).tolist()

        output = [score_sum / pos_sum for score_sum, pos_sum in zip(scores_summed, pos_summed)]
        return output


def main():
    args = arg_parse()
    avgCalc = AvgCalc(args)
    if args.chunker:
        with sys.stdin as chunk:
            x, y = avgCalc.chunk_handler(chunk)
            print(list(x))
            print(list(y))

    elif args.combine:
        scores = []
        pos = []
        with sys.stdin as res:
            for i, line in enumerate(res.readlines()):
                if (i + 1) % 2 != 0:
                    scores.append(list(map(float, line.strip().strip("[] ").replace(" ", "").split(","))))
                else:
                    pos.append(list(map(float, line.strip().strip("[] ").replace(" ", "").split(","))))

        scores = [np.array(result) for result in scores]
        pos = [np.array(result) for result in pos]

        max_len = len(max(scores, key=len))

        arr1 = np.array([np.pad(row, (0, max_len - len(row)), 'constant') for row in scores])
        arr2 = np.array([np.pad(row, (0, max_len - len(row)), 'constant') for row in pos])

        scores_summed = np.sum(arr1, axis=0).tolist()
        pos_summed = np.sum(arr2, axis=0).tolist()

        output = [score_sum / pos_sum for score_sum, pos_sum in zip(scores_summed, pos_summed)]

        df = pd.DataFrame(output)
        df.to_csv(sys.stdout, header=False)

if __name__ == "__main__":
    main()
