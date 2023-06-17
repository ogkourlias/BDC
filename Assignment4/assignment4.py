#!/usr/bin/env python3

import argparse as ap
import sys
from pathlib import Path
from mpi4py import MPI
import pandas as pd
import numpy as np

comm = MPI.COMM_WORLD
comm_size = comm.Get_size()
my_rank = comm.Get_rank()

def parse_arguments():
    argparser = ap.ArgumentParser(description="Calculate PHRED scores over the network.")
    server_args = argparser.add_argument_group(title="Arguments for server mode")
    server_args.add_argument("-o", dest="csvfile", type=Path, help="CSV file to store the output. Defaults to STDOUT.")
    server_args.add_argument("fastq_files", type=Path, nargs='*', help="At least one FastQ file to process.")
    return argparser.parse_args()

class PhredScoreCalculator:

    def get_score_lines(self, fastq):
        lines = fastq.readlines()
        return [line.strip() for i, line in enumerate(lines) if (i + 1) % 4 == 0]

    def divide_work(self, score_lines, chunk_size):
        chunk_size = int(chunk_size)
        return [score_lines[i:i + chunk_size] for i in range(0, len(score_lines), chunk_size)]

    def get_scores(self, score_lines):
        scores = []
        counts = []
        for line in score_lines:
            for i, score in enumerate(line):
                score_value = ord(score) - 33
                if i < len(scores):
                    scores[i] += score_value
                    counts[i] += 1
                else:
                    scores.append(score_value)
                    counts.append(1)
        return scores, counts

    def combine_results(self, results):
        max_len = max(len(result[0]) for result in results)
        scores_arrays = [np.pad(result[0], (0, max_len - len(result[0])), 'constant') for result in results]
        counts_arrays = [np.pad(result[1], (0, max_len - len(result[1])), 'constant') for result in results]
        scores = np.sum(scores_arrays, axis=0).tolist()
        counts = np.sum(counts_arrays, axis=0).tolist()
        return [score / count for score, count in zip(scores, counts)]

    def save_to_csv(self, scores, csvfile = None):
        df = pd.DataFrame(scores, columns=['Scores'])
        if csvfile:
            df.to_csv(csvfile, header=False)
        else:
            df.to_csv(sys.stdout, header=False)

    def files_handler(self, files):
        file_line_list = []
        for file in files:
            with open(file) as fastq:
                lines = fastq.readlines()
                file_line_list.append([line.strip() for i, line in enumerate(lines) if (i + 1) % 4 == 0])
        return file_line_list

def main():
    args = parse_arguments()
    calculator = PhredScoreCalculator()
    if my_rank == 0:
        data = calculator.files_handler(args.fastq_files)
        work = calculator.divide_work(data[0], (len(data[0]) / comm_size))
    else:
        work = None
    work_this = comm.scatter(work, root=0)
    scores = [calculator.get_scores(work_this)]
    res = calculator.combine_results(scores)
    response_array = comm.gather(res, root=0)
    if my_rank == 0:
        calculator.save_to_csv(response_array[0])

if __name__ == "__main__":
    main()
