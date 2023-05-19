# Heb het eerst met een class implementatie geprobeerd. Ik kreeg hier helaas de hele tijd een
# error. Volgens mij omdat pool top-level functies nodig heeft.
class QualityController:
    """ Desc """

    def __init__(self):
        # self.fastqfile = '/commons/Themas/Thema12/HPC/rnaseq_selection.fastq'
        self.args = self.arg_parse()
        self.linelength = 101
        self.scores = [0 for _ in range(self.linelength)]
        self.cpus = self.args.n
        self.files = self.args.fastq_files
        print(self.files)
        self.read_count = 0

    def arg_parse(self):
        argparser = ap.ArgumentParser(description="Script voor Opdracht 1 van Big Data Computing")
        argparser.add_argument("-n", action="store",
                               dest="n", required=True, type=int,
                               help="Aantal cores om te gebruiken.")
        argparser.add_argument("-o", action="store", dest="csvfile", type=ap.FileType('w', encoding='UTF-8'),
                               required=False,
                               help="CSV file om de output in op te slaan. Default is output naar terminal STDOUT")
        argparser.add_argument("fastq_files", action="store", nargs='+',
                               help="Minstens 1 Illumina Fastq Format file om te verwerken")
        return argparser.parse_args()

    def avg_accuracy(self, fastq):
        """ Iterate through fastqfile. NB: chunk is assumed to be multiple of 4!"""
        read_count = 0
        with open(fastq, 'r') as fastq:
            for i, line in enumerate(fastq):
                if (i + 1) % 4 == 0:
                    read_count += 1
                    for score_i in range(len(self.scores)):
                        self.scores[score_i] += 1 - 10 ** -((ord(line[score_i]) - 33) / 10)

        print()

        return [score / read_count for score in self.scores]

    def handler(self):
        with Pool(self.cpus) as p:
            res = p.map(self.avg_accuracy, self.files)