import multiprocessing as mp
import argparse as ap
import queue
import sys
from demo import *
import pandas as pd

def argparser():
    """Parse command line arguments."""
    argparser = ap.ArgumentParser(
        description="Script voor Opdracht 2 van Big Data Computing;  Calculate PHRED scores over the network.")
    mode = argparser.add_mutually_exclusive_group(required=True)
    mode.add_argument("-s", action="store_true", help="Run the program in Server mode; see extra options needed below")
    mode.add_argument("-c", action="store_true", help="Run the program in Client mode; see extra options needed below")
    server_args = argparser.add_argument_group(title="Arguments when run in server mode")
    server_args.add_argument("-o", action="store", dest="csvfile", type=ap.FileType('w', encoding='UTF-8'),
                             required=False,
                             help="CSV file om de output in op te slaan. Default is output naar terminal STDOUT")
    server_args.add_argument("fastq_files", action="store", type=ap.FileType('r'), nargs='*',
                             help="Minstens 1 Illumina Fastq Format file om te verwerken")
    server_args.add_argument("--chunks", action="store", type=int, required=True)

    client_args = argparser.add_argument_group(title="Arguments when run in client mode")
    client_args.add_argument("-n", action="store",
                             dest="n", required=False, type=int,
                             help="Aantal cores om te gebruiken per host.")
    client_args.add_argument("--host", action="store", type=str, help="The hostname where the Server is listening")
    client_args.add_argument("--port", action="store", type=int, help="The port on which the Server is listening")

    args = argparser.parse_args()
    return args

def server(file, chunksize = 50000):
    manager = make_server_manager(AUTHKEY, PORTNUM)
    q_man = manager.join()
    read_count = 0
    with file as fastq:
        lines = fastq.readlines()
        score_lines = []
        for i, line in enumerate(lines):
            if (i + 1) % 4 == 0:
                score_lines.append(line.strip())
                read_count += 1
                jobs = [score_lines[i:i+chunksize] for i in range(0, len(score_lines), chunksize)]
    for job in jobs:
        man

    manager.


def walk(file, cores):
    """Iterate through fastqfile. NB: chunk is assumed to be multiple of 4!"""
    pool = mp.Pool(cores)
    read_count = 0
    chunksize = 500_000
    with open(file) as fastq:
        lines = fastq.readlines()
        score_lines = []
        for i, line in enumerate(lines):
            if (i + 1) % 4 == 0:
                score_lines.append(line.strip())
                read_count += 1

        scores = pool.map(line_handler, [score_lines[i:i+chunksize]
                                         for i in range(0, len(score_lines), chunksize)])
        scores = [sum(score)/read_count for score in zip(*scores)]

        return scores

def line_handler(lines):
    """Calculate average accuracy for a chunk of lines."""
    scores = [0 for _ in range(101)]
    for line in lines:
        for score_i in range(len(scores)):
            scores[score_i] += (ord(line[score_i]) - 33)
    return scores

def write_csv(scores, csv_file):
    """Write scores to csv file."""
    df = pd.DataFrame(scores)
    df.to_csv(csv_file, header=False)

def write_to_stdout(scores):
    """Write scores to stdout."""
    df = pd.DataFrame(scores)
    df.to_csv(sys.stdout, header=False)

def calc_avg(scores, read_count):
    """Calculate average accuracy."""
    return [score / read_count for score in scores]

def main():
    """Main function."""
    args = argparser()
    scores = walk(args.fastq_files[0], args.n)
    if args.csvfile:
        write_csv(scores, args.csvfile)
    else:
        write_to_stdout(scores)


if __name__ == "__main__":
    main()