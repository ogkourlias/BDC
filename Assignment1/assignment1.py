#!/usr/bin/env python3

"""python3 assignment1.py -n <aantal_cpus> [OPTIONEEL: -o <output csv file>] fastabestand1.fastq
[fastabestand2.fastq ... fastabestandN.fastq]"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.1"

# IMPORTS
import argparse as ap
import csv
import sys
from multiprocessing import Pool


# FUNCTIONS
def arg_parse():
    """Parse arguments"""
    argparser = ap.ArgumentParser(
        description="Script voor Opdracht 1 van Big Data Computing"
    )
    argparser.add_argument(
        "-n",
        action="store",
        dest="n",
        required=True,
        type=int,
        help="Aantal cores om te gebruiken.",
    )
    argparser.add_argument(
        "-o",
        action="store",
        dest="csvfile",
        required=False,
        help="CSV file om de output in op te slaan. Default is output naar terminal STDOUT",
    )
    argparser.add_argument(
        "fastq_files",
        action="store",
        nargs="+",
        help="Minstens 1 Illumina Fastq Format file om te verwerken",
    )
    return argparser.parse_args()


def quality_avg_jobs(cpus, files):
    """Create jobs for multiprocessing."""
    pool = Pool(cpus)
    res = pool.map(avg_accuracy, files)
    pool.close()
    pool.join()
    return res


def avg_accuracy(fastq):
    """Iterate through fastqfile. NB: chunk is assumed to be multiple of 4!"""
    linelength = 101
    scores = [0 for _ in range(linelength)]
    read_count = 0
    with open(fastq, "r") as fastq:
        for i, line in enumerate(fastq):
            if (i + 1) % 4 == 0:
                read_count += 1
                for score_i in range(len(scores)):
                    scores[score_i] += 1 - 10 ** -((ord(line[score_i]) - 33) / 10)

    return [score / read_count for score in scores]


def output_jobs(cpus, results, o, files):
    """Create output jobs to be executed."""
    if o:
        for file, result in zip(files, results):
            write_csv(result, o, file)
    else:
        for file, result in zip(files, results):
            sys.stdout.write(f"{file.split('/')[-1]}\n{result}\n")


def write_csv(result, output, file):
    """Write results to csv file."""
    with open(output, "a+") as csvfile:
        writer = csv.writer(csvfile)
        csvfile.write(file.split("/")[-1] + "\n")
        writer.writerow(result)


# MAIN
def main():
    """Main function"""
    # Argparse
    args = arg_parse()
    cpus = args.n
    files = args.fastq_files
    output = args.csvfile

    # Multiprocessing function calls
    res = quality_avg_jobs(cpus, files)
    output_jobs(cpus, res, output, files)
    # FINISH
    return 0


if __name__ == "__main__":
    main()
