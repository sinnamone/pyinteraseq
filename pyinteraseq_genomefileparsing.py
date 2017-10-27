import os
from Bio import SeqIO
import pandas as pd
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
            self.chromosomename = self.seq_record.id.split("|")[-2]  # TODO richiede stringa non modificare file fasta
            self.chromosomelength = len(self.seq_record)

    def genomefilewrite(self):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            with open(self.outputfolder + self.chromosomename + '.genome', "wb") as self.genome:
                self.genome.write(self.chromosomename + '\t' + str(self.chromosomelength))
        except OSError:
            self.filelog.write(msg45)
            sys.exit(0)
        else:
            self.filelog.write(msg46)
            return self.outputfolder + self.chromosomename + '.genome'

    def fastareference(self):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            with open(self.outputfolder + self.chromosomename + '.fasta', 'w') as self.ref:
                for self.seq_record in SeqIO.parse(self.fastasequence, "fasta"):
                    self.ref.write('>' + str(self.chromosomename) + '\n' + str(self.seq_record.seq))
        except OSError:
            self.filelog.write(msg43)
            sys.exit(0)
        else:
            for item in self.dir:
                if item.endswith(".fai"):
                    os.remove(self.outputfolder+item)
            self.filelog.write(msg44)
            self.index = pysam.Fastafile(self.outputfolder + self.chromosomename + '.fasta')
            return self.outputfolder + self.chromosomename + '.fasta'


class AnnotationFile(GenomeFile):

    def __init__(self, optparseinstance):
        GenomeFile.__init__(self, optparseinstance)
        self.NameFileAnnotation = os.path.basename(self.annotation.split('/')[-1])
        self.PathAnnotation = os.path.dirname(self.annotation)
        self.dfAnno = None
        self.aux = None

    def annotationbuild(self):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            self.dfAnno = pd.read_csv(self.annotation, sep="\t", header=2)
            self.aux = self.dfAnno['Location'].apply(lambda x: x.split('..'))
            self.dfAnno['Start'] = self.aux.apply(lambda x: x[0])
            self.dfAnno['End'] = self.aux.apply(lambda x: x[1])
            self.dfAnno['Chr'] = self.chromosomename
            self.dfAnno[['Chr', 'Start', 'End', 'Synonym', 'COG', 'Strand', 'Gene',
                         'Product']].to_csv(self.outputfolder + self.chromosomename + '_proteome.bed',
                                            header=None, sep='\t', index=False)
        except OSError:
            self.filelog.write(msg41)
            sys.exit(0)
        else:
            self.filelog.write(msg42)
            return self.outputfolder + self.chromosomename + '_proteome.bed'
