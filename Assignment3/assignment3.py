#!/usr/bin/env python3

"""
As Server:
assignment2.py -s -ch <chunksize> -i <input_fasta_path> -o <output_csv_path>

As Client:
assignment2.py -c -n <amount_of_cores>"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.1"

# IMPORTS
import sys, os
import argparse as ap
from pathlib import Path
import multiprocessing as mp
import pandas as pd
import numpy as np
from itertools import zip_longest


# FUNCTIONS
def arg_parse():
    # Generic args.
    argparser = ap.ArgumentParser(
        description="Script voor Opdracht 2 van Big Data Computing;  Calculate PHRED scores over the network."
    )
    mode = argparser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "-s",
        action="store_true",
        help="Run the program in Server mode; see extra options needed below",
    )
    mode.add_argument(
        "-c",
        action="store_true",
        help="Run the program in Client mode; see extra options needed below",
    )
    mode.add_argument(
        "-calc",
        action="store_true",
        help="Run the program in calc mode; see extra options needed below",
    )

    # Server args.
    server_args = argparser.add_argument_group(
        title="Arguments when run in server mode"
    )

    server_args.add_argument(
        "-o",
        action="store",
        dest="output",
        type=Path,
        required=False,
        help="CSV file om de output in op te slaan. Default is output naar terminal STDOUT",
    )

    server_args.add_argument(
        "-i",
        action="store",
        dest="input",
        type=Path,
        nargs="*",
        help="Minstens 1 Illumina Fastq Format file om te verwerken",
    )

    server_args.add_argument("-ch", action="store", dest="chunksize", type=int)

    # Client args.
    client_args = argparser.add_argument_group(
        title="Arguments when run in client mode"
    )
    client_args.add_argument(
        "-n",
        action="store",
        dest="n",
        required=False,
        type=int,
        help="Aantal cores om te gebruiken per host.",
    )

    client_args.add_argument(
        "-seq",
        action="store",
        dest="seq",
        required=False,
        type=Path,
        help="Aantal cores om te gebruiken per host.",
    )

    calc_args = argparser.add_argument_group(title="Arguments when run in client mode")
    calc_args.add_argument(
        "-co",
        action="store",
        dest="calc_out",
        required=False,
        type=str,
        help="Aantal cores om te gebruiken per host.",
    )

    calc_args.add_argument(
        "-if",
        action="store",
        dest="calc_input",
        required=False,
        type=list,
        help="Aantal cores om te gebruiken per host.",
    )

    # Calc args.
    return argparser.parse_args()


class AvgCalc:
    """
    Class containg functiosn to calculate the average Phred scores for positions.
    """

    def __init__(self, args):
        self.files = args.input
        self.chunk_size = args.chunksize
        self.cores = args.n
        self.fastq_files = args.input
        self.csvfile = args.output

    def calculate(self, files):
        """
        Function containing the primary calculation pipeline for  given files.
        """
        for file in files:
            with file as fastq:
                self.line_walker(fastq)

    def line_walker(self, file, cores):
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

            res = pool.map(
                self.score_getter,
                [
                    score_lines[i : i + self.chunk_size]
                    for i in range(0, len(score_lines), self.chunk_size)
                ],
            )

            score_list = [np.array(result[0]) for result in res]
            pos_list = [np.array(result[1]) for result in res]
            arr1 = np.array(
                [
                    np.pad(
                        row, (0, len(max(score_list, key=len)) - len(row)), "constant"
                    )
                    for row in score_list
                ]
            )
            arr2 = np.array(
                [
                    np.pad(
                        row, (0, len(max(score_list, key=len)) - len(row)), "constant"
                    )
                    for row in pos_list
                ]
            )
            scores_summed = np.sum(arr1, axis=0).tolist()
            pos_summed = np.sum(arr2, axis=0).tolist()
            output = [
                score_sum / pos_sum
                for score_sum, pos_sum in zip(scores_summed, pos_summed)
            ]

    def line_walker_client(lines):
        """
        Walks through the lines the given input
        :param file:
        """
        pos_scores = []
        pos_counts = []
        pos_dict = {}
        for line in lines:
            for i in range(len(line)):
                if len(pos_scores) > i:
                    pos_dict[i] = pos_dict[i] + 1
                    pos_counts[i] += 1
                    pos_scores[i] += ord(line[i]) - 33
                else:
                    pos_dict[i] = 1
                    pos_counts.append(1)
                    pos_scores.append(ord(line[i]) - 33)

        return pos_scores

    def get_res_server(res):
        score_list = [np.array(result[0]) for result in res]
        pos_list = [np.array(result[1]) for result in res]
        arr1 = np.array(
            [
                np.pad(row, (0, len(max(score_list, key=len)) - len(row)), "constant")
                for row in score_list
            ]
        )
        arr2 = np.array(
            [
                np.pad(row, (0, len(max(score_list, key=len)) - len(row)), "constant")
                for row in pos_list
            ]
        )
        scores_summed = np.sum(arr1, axis=0).tolist()
        pos_summed = np.sum(arr2, axis=0).tolist()
        output = [
            score_sum / pos_sum for score_sum, pos_sum in zip(scores_summed, pos_summed)
        ]
        return output

    def read_lines(file):
        """
        Walks through the lines the given input
        :param file:
        """
        with open(file) as fastq:
            score_lines = []
            lines = fastq.readlines()
            for i, line in enumerate(lines):
                if (i + 1) % 4 == 0:
                    score_lines.append(line.strip())
            return score_lines

    def work_division(self, lines, chunk_size):
        res = [lines[i : i + chunk_size] for i in range(0, len(lines), chunk_size)]
        return res

    def work_division_server(lines, chunk_size):
        res = [lines[i : i + chunk_size] for i in range(0, len(lines), chunk_size)]
        return res

    def calc_per_file(self, work):
        lines_for_each = work
        cores = self.cores
        chunk = self.chunk_size
        pool = mp.Pool(cores)
        final = []
        for line_list in lines_for_each:
            work = self.work_division(line_list, chunk)
            res = pool.map(self.score_getter, work)
            score_list = [np.array(result[0]) for result in res]
            pos_list = [np.array(result[1]) for result in res]
            arr1 = np.array(
                [
                    np.pad(
                        row, (0, len(max(score_list, key=len)) - len(row)), "constant"
                    )
                    for row in score_list
                ]
            )
            arr2 = np.array(
                [
                    np.pad(
                        row, (0, len(max(score_list, key=len)) - len(row)), "constant"
                    )
                    for row in pos_list
                ]
            )
            scores_summed = np.sum(arr1, axis=0).tolist()
            pos_summed = np.sum(arr2, axis=0).tolist()
            output = [
                score_sum / pos_sum
                for score_sum, pos_sum in zip(scores_summed, pos_summed)
            ]
            final.append(output)

        return final

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
                    pos_scores[i] += ord(line[i]) - 33
                else:
                    pos_dict[i] = 1
                    pos_counts.append(1)
                    pos_scores.append(ord(line[i]) - 33)

        return pos_scores, pos_counts

    def score_getter_client(lines):
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
                    pos_scores[i] += ord(line[i]) - 33
                else:
                    pos_dict[i] = 1
                    pos_counts.append(1)
                    pos_scores.append(ord(line[i]) - 33)

        return pos_scores, pos_counts

    def score_getter_line(line):
        """
        Gets the scores and indexes for lines
        :param lines:
        :return:
        """
        pos_scores = []
        pos_counts = []
        pos_dict = {}
        for i in range(len(line)):
            if len(pos_scores) > i:
                pos_dict[i] = pos_dict[i] + 1
                pos_counts[i] += 1
                pos_scores[i] += ord(line[i]) - 33
            else:
                pos_dict[i] = 1
                pos_counts.append(1)
                pos_scores.append(ord(line[i]) - 33)

        return pos_scores

    def line_worker(self, work):
        res = self.score_getter(work)
        score_list = [np.array(result[0]) for result in res]
        pos_list = [np.array(result[1]) for result in res]
        arr1 = np.array(
            [
                np.pad(row, (0, len(max(score_list, key=len)) - len(row)), "constant")
                for row in score_list
            ]
        )
        arr2 = np.array(
            [
                np.pad(row, (0, len(max(score_list, key=len)) - len(row)), "constant")
                for row in pos_list
            ]
        )
        scores_summed = np.sum(arr1, axis=0).tolist()
        pos_summed = np.sum(arr2, axis=0).tolist()
        output = [
            score_sum / pos_sum for score_sum, pos_sum in zip(scores_summed, pos_summed)
        ]
        return output

    def files_handler(self, files):
        file_line_list = []
        for file in files:
            with open(file) as fastq:
                lines = fastq.readlines()
                file_line_list.append(
                    [line.strip() for i, line in enumerate(lines) if (i + 1) % 4 == 0]
                )

        return file_line_list

    def calc(self, res):
        score_list = [np.array(result[0]) for result in res]
        pos_list = [np.array(result[1]) for result in res]
        arr1 = np.array(
            [
                np.pad(row, (0, len(max(score_list, key=len)) - len(row)), "constant")
                for row in score_list
            ]
        )
        arr2 = np.array(
            [
                np.pad(row, (0, len(max(score_list, key=len)) - len(row)), "constant")
                for row in pos_list
            ]
        )
        scores_summed = np.sum(arr1, axis=0).tolist()
        pos_summed = np.sum(arr2, axis=0).tolist()
        output = [
            score_sum / pos_sum for score_sum, pos_sum in zip(scores_summed, pos_summed)
        ]
        return output

    def calc_from_file(score_list, pos_list):
        print(score_list)
        arr1 = np.array(
            [
                np.pad(row, (0, len(max(score_list, key=len)) - len(row)), "constant")
                for row in score_list
            ]
        )
        arr2 = np.array(
            [
                np.pad(row, (0, len(max(score_list, key=len)) - len(row)), "constant")
                for row in pos_list
            ]
        )
        scores_summed = np.sum(arr1, axis=0).tolist()
        pos_summed = np.sum(arr2, axis=0).tolist()
        output = [
            score_sum / pos_sum for score_sum, pos_sum in zip(scores_summed, pos_summed)
        ]
        return output

    def write_output(self, output):
        """Write scores to csv file."""
        df = pd.DataFrame(output)
        if self.csvfile:
            df.to_csv(self.csvfile, header=False)
        else:
            df.to_csv(sys.stdout, header=False)

    def write_output_server(output, output_file):
        """Write scores to csv file."""
        df = pd.DataFrame(output)
        if output_file:
            df.to_csv(output_file, header=False)
        else:
            df.to_csv(sys.stdout, header=False)


# MAIN
def main():
    """Main function"""
    args = arg_parse()
    if args.c:
        with open(args.seq, "r+") as f:
            with open(f"scores/{str(args.seq)[4:]}_output", "a+") as output_f:
                for line in f:
                    output_f.write(f"{AvgCalc.score_getter_line(line.strip())}\n")

    elif args.s:
        lines = []
        counter = 0
        with open(args.input[0]) as f:
            for i, line in enumerate(f):
                if (i + 1) % 4 == 0:
                    lines.append(line)
                if len(lines) == args.chunksize:
                    with open(f"work/work_{counter}.txt", "w+") as work_f:
                        for line in lines:
                            work_f.write(f"{line}")
                    lines = []
                    counter += 1

        if len(lines) != 0:
            with open(f"work/work_{counter}", "w+") as work_f:
                for line in lines:
                    work_f.write(f"{line}")

    elif args.calc:
        files = os.listdir("scores/")
        files = [f"scores/{file}" for file in files]
        sums = []
        positions = []
        for file in files:
            with open(file, "r+") as f:
                for line in f:
                    line = line[1:-2].split(",")
                    scores = [int(score) for score in line]
                    for i, (score, summed, pos) in enumerate(
                        zip_longest(scores, sums, positions, fillvalue=0)
                    ):
                        if i >= len(sums):
                            sums.append(summed + score)
                            positions.append(pos + 1)
                        else:
                            sums[i] += score
                            if score != 0:
                                positions[i] += 1

        output = [score_sum / pos_sum for score_sum, pos_sum in zip(sums, positions)]
        AvgCalc.write_output_server(output, args.calc_out)


if __name__ == "__main__":
    main()
