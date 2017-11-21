import os
from Bio import SeqIO
import sys
from output_message import *
import pysam
from pyinteraseq_inputcheck import InputCheck


class GenomeFile(InputCheck):

    def __init__(self, optparseinstance):
        InputCheck.__init__(self, optparseinstance)
        self.namefilefasta = os.path.basename(self.fastasequence.split('/')[-1])
        self.PathFasta = os.path.dirname(self.fastasequence)
        self.genome = None
        self.ref = None
        self.index = None
        self.dir = os.listdir(self.outputfolder)
        for self.seq_record in SeqIO.parse(self.fastasequence, "fasta"):
            self.chromosomelength = len(self.seq_record)

    def genomefilewrite(self):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            with open(self.outputfolder + self.genename + '.genome', "wb") as self.genome:
                self.genome.write(self.genename + '\t' + str(self.chromosomelength))
        except OSError:
            self.filelog.write(msg45)
            sys.exit(0)
        else:
            self.filelog.write(msg46)
            return self.outputfolder + self.genename + '.genome'

    def fastareference(self):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            with open(self.outputfolder + self.genename + '.fasta', 'w') as self.ref:
                for self.seq_record in SeqIO.parse(self.fastasequence, "fasta"):
                    self.ref.write('>' + str(self.genename) + '\n' + str(self.seq_record.seq))
        except OSError:
            self.filelog.write(msg43)
            sys.exit(0)
        else:
            for item in self.dir:
                if item.endswith(".fai"):
                    os.remove(self.outputfolder+item)
            self.filelog.write(msg44)
            self.index = pysam.Fastafile(self.outputfolder + self.genename + '.fasta')
            return self.outputfolder + self.genename + '.fasta'







