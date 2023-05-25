import multiprocessing as mp
import argparse as ap
import sys

import pandas as pd

def argparser():
    argparser = ap.ArgumentParser(description="Script voor Opdracht 1 van Big Data Computing")
    argparser.add_argument("-n", action="store",
                           dest="n", required=True, type=int,
                           help="Aantal cores om te gebruiken.")
    argparser.add_argument("-o", action="store", dest="csvfile", type=ap.FileType('w', encoding='UTF-8'),
                           required=False,
                           help="CSV file om de output in op te slaan. Default is output naar terminal STDOUT")
    argparser.add_argument("fastq_files", action="store", nargs='+',
                           help="Minstens 1 Illumina Fastq Format file om te verwerken")
    args = argparser.parse_args()
    return args

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

        scores = pool.map(line_handler, [score_lines[i:i+chunksize] for i in range(0, len(score_lines), chunksize)])
        scores = [sum(score)/read_count for score in zip(*scores)]

        return scores

def line_handler(lines):
    scores = [0 for _ in range(101)]
    for line in lines:
        for score_i in range(len(scores)):
            scores[score_i] += (ord(line[score_i]) - 33)
    return scores

def write_csv(scores, csv_file):
    df = pd.DataFrame(scores)
    df.to_csv(csv_file, header=False)

def write_to_stdout(scores):
    df = pd.DataFrame(scores)
    df.to_csv(sys.stdout, header=False)

def calc_avg(scores, read_count):
    return [score / read_count for score in scores]

def main():
    args = argparser()
    scores = walk(args.fastq_files[0], args.n)
    if args.csvfile:
        write_csv(scores, args.csvfile)
    else:
        write_to_stdout(scores)

if __name__ == "__main__":
    main()