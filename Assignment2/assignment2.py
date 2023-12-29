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
POISONPILL = "MEMENTOMORI"
ERROR = "DOH"
IP = ''
PORTNUM = 5381
AUTHKEY = b'whathasitgotinitspocketsesss?'
data = ["Always", "look", "on", "the", "bright", "side", "of", "life!"]


def make_server_manager(port, authkey):
    """ Create a manager for the server, listening on the given port.
        Return a manager object with get_job_q and get_result_q methods.
    """
    job_q = queue.Queue()
    result_q = queue.Queue()

    # This is based on the examples in the official docs of multiprocessing.
    # get_{job|result}_q return synchronized proxies for the actual Queue
    # objects.
    class QueueManager(BaseManager):
        pass

    QueueManager.register('get_job_q', callable=lambda: job_q)
    QueueManager.register('get_result_q', callable=lambda: result_q)

    manager = QueueManager(address=('', port), authkey=authkey)
    manager.start()
    print('Server started at port %s' % port)
    return manager


def runserver(fn, data):
    # Start a shared manager server and access its queues
    manager = make_server_manager(PORTNUM, b'whathasitgotinitspocketsesss?')
    shared_job_q = manager.get_job_q()
    shared_result_q = manager.get_result_q()

    if not data:
        print("Gimme something to do here!")
        return
    print("Sending data!")
    for d in data:
        shared_job_q.put({'fn': fn, 'arg': d})
    time.sleep(2)
    results = []
    while True:
        try:
            result = shared_result_q.get_nowait()
            results.append(result)
            print("Got result!", result)
            if len(results) == len(data):
                print("Got all results!")
                break
        except queue.Empty:
            time.sleep(1)
            continue
    # Tell the client process no more data will be forthcoming
    print("Time to kill some peons!")
    shared_job_q.put(POISONPILL)
    # Sleep a bit before shutting down the server - to give clients time to
    # realize the job queue is empty and exit in an orderly way.
    time.sleep(5)
    print("Aaaaaand we're done for the server!")
    manager.shutdown()
    print(results)


def make_client_manager(ip, port, authkey):
    """ Create a manager for a client. This manager connects to a server on the
        given address and exposes the get_job_q and get_result_q methods for
        accessing the shared queues from the server.
        Return a manager object.
    """
    class ServerQueueManager(BaseManager):
        pass

    ServerQueueManager.register('get_job_q')
    ServerQueueManager.register('get_result_q')

    manager = ServerQueueManager(address=(ip, port), authkey=authkey)
    manager.connect()

    print('Client connected to %s:%s' % (ip, port))
    return manager


def capitalize(word):
    """Capitalizes the word you pass in and returns it"""
    return word.upper()


def runclient(num_processes):
    manager = make_client_manager(IP, PORTNUM, AUTHKEY)
    job_q = manager.get_job_q()
    result_q = manager.get_result_q()
    run_workers(job_q, result_q, num_processes)


def run_workers(job_q, result_q, num_processes):
    processes = []
    for p in range(num_processes):
        temP = mp.Process(target=peon, args=(job_q, result_q))
        processes.append(temP)
        temP.start()
    print("Started %s workers!" % len(processes))
    for temP in processes:
        temP.join()


def peon(job_q, result_q):
    my_name = mp.current_process().name
    while True:
        try:
            job = job_q.get_nowait()
            if job == POISONPILL:
                job_q.put(POISONPILL)
                print("Aaaaaaargh", my_name)
                return
            else:
                try:
                    result = job['fn'](job['arg'])
                    print("Peon %s Workwork on %s!" % (my_name, job['arg']))
                    result_q.put({'job': job, 'result': result})
                except NameError:
                    print("Can't find yer fun Bob!")
                    result_q.put({'job': job, 'result': ERROR})

        except queue.Empty:
            print("sleepytime for", my_name)
            time.sleep(1)


def arg_parse():
    # Generic args.
    argparser = ap.ArgumentParser(
        description="Script voor Opdracht 2 van Big Data Computing;  Calculate PHRED scores over the network.")
    mode = argparser.add_mutually_exclusive_group(required=True)
    mode.add_argument("-s", action="store_true", help="Run the program in Server mode; see extra options needed below")
    mode.add_argument("-c", action="store_true", help="Run the program in Client mode; see extra options needed below")
    
    # Server args.
    server_args = argparser.add_argument_group(title="Arguments when run in server mode")
    
    server_args.add_argument("-o", action="store", dest="output", type=Path,
                             required=False,
                             help="CSV file om de output in op te slaan. Default is output naar terminal STDOUT")
    
    server_args.add_argument("-i", action="store", dest="input", type=Path, nargs='*',
                             help="Minstens 1 Illumina Fastq Format file om te verwerken")
    
    server_args.add_argument("-ch", action="store", dest="chunksize", type=int, required=True)

    # Client args.
    client_args = argparser.add_argument_group(title="Arguments when run in client mode")
    client_args.add_argument("-n", action="store",
                             dest="n", required=False, type=int,
                             help="Aantal cores om te gebruiken per host.")
    
    client_args.add_argument("-hn", action="store", dest="host", type=str, help="The hostname where the Server is listening")
    
    client_args.add_argument("-p", action="store", dest="port",type=int, help="The port on which the Server is listening")
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

    def work_division(self, lines, chunk_size):
        res = [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)]
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
                [np.pad(row, (0, len(max(score_list, key=len)) - len(row)), 'constant') for row in score_list])
            arr2 = np.array(
                [np.pad(row, (0, len(max(score_list, key=len)) - len(row)), 'constant') for row in pos_list])
            scores_summed = np.sum(arr1, axis=0).tolist()
            pos_summed = np.sum(arr2, axis=0).tolist()
            output = [score_sum / pos_sum for score_sum, pos_sum in zip(scores_summed, pos_summed)]
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
                    pos_scores[i] += (ord(line[i]) - 33)
                else:
                    pos_dict[i] = 1
                    pos_counts.append(1)
                    pos_scores.append(ord(line[i]) - 33)

        return pos_scores, pos_counts

    def line_worker(self, work):
        res = self.score_getter(work)
        score_list = [np.array(result[0]) for result in res]
        pos_list = [np.array(result[1]) for result in res]
        arr1 = np.array(
            [np.pad(row, (0, len(max(score_list, key=len)) - len(row)), 'constant') for row in score_list])
        arr2 = np.array(
            [np.pad(row, (0, len(max(score_list, key=len)) - len(row)), 'constant') for row in pos_list])
        scores_summed = np.sum(arr1, axis=0).tolist()
        pos_summed = np.sum(arr2, axis=0).tolist()
        output = [score_sum / pos_sum for score_sum, pos_sum in zip(scores_summed, pos_summed)]
        return output

    def files_handler(self, files):
        file_line_list = []
        for file in files:
            with open(file) as fastq:
                lines = fastq.readlines()
                file_line_list.append([line.strip() for i, line in enumerate(lines) if (i + 1) % 4 == 0])

        return file_line_list

    def calc(self, res):
        score_list = [np.array(result[0]) for result in res]
        pos_list = [np.array(result[1]) for result in res]
        arr1 = np.array(
            [np.pad(row, (0, len(max(score_list, key=len)) - len(row)), 'constant') for row in score_list])
        arr2 = np.array(
            [np.pad(row, (0, len(max(score_list, key=len)) - len(row)), 'constant') for row in pos_list])
        scores_summed = np.sum(arr1, axis=0).tolist()
        pos_summed = np.sum(arr2, axis=0).tolist()
        output = [score_sum / pos_sum for score_sum, pos_sum in zip(scores_summed, pos_summed)]
        return output

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
    avgCalc = AvgCalc(args)
    if args.c:
        client = mp.Process(target=runclient, args=(4,))
        client.start()
        server.join()
        client.join()
    elif args.s:
        server = mp.Process(target=runserver, args=(capitalize, data))
        server.start()
        time.sleep(1)
        print(args)
    
    # server = mp.Process(target=runserver, args=(avgCalc.calc_per_file, work))
    # server.start()
    # time.sleep(1)
    # client = mp.Process(target=runclient, args=(args.n,))
    # client.start()
    # server.join()
    # client.join()

if __name__ == "__main__":
    main()
