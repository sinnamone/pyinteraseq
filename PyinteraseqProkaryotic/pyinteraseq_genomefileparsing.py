import os
from Bio import SeqIO
import pandas as pd
import sys
from output_message import *
import pysam
import pybedtools
# import warnings
import traceback
from pyinteraseq_inputcheck import InputCheck


class GenomeFile(InputCheck):

    def __init__(self, optparseinstance):
        InputCheck.__init__(self, optparseinstance)
        #name of fasta file without path
        self.namefilefasta = os.path.basename(self.fastasequence.split('/')[-1])
        # path were fasta file is located
        self.PathFasta = os.path.dirname(self.fastasequence)
        self.chromosomename = self.namefilefasta.split('.')[0]
        self.genome = None
        self.ref = None
        self.index = None
        self.numchr = None
        self.df1 = None
        self.dir = os.listdir(self.outputfolder)
        self.id = '1'
        self.id = [(int(x) - 1) for x in self.id]
        self.header = None
        self.seqix = 1
        for self.seq_record in SeqIO.parse(self.fastasequence, "fasta"):
            # save the name of chromosome. Format string is composed by >gi|XXXXXX|ref|NC_XXXXX.1|description
            # self.chromosomename = self.seq_record.id.split("|")[-2]
            self.chromosomelength = len(self.seq_record)
        print self.chromosomename


    def genomefilewrite(self):
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
        try:
            with open(self.outputfolder + self.chromosomename + '.genome', "wb") as self.genome:
                print self.genome
                self.genome.write(self.chromosomename + '\t' + str(self.chromosomelength))
        except OSError:
            self.filelog.write(msg45)
            sys.exit(0)
        else:
            self.filelog.write(msg46)
            return self.outputfolder + self.chromosomename + '.genome'

    def ref2tabular(self):
        """
        Function for conversion Reference Fasta file in tabular format
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
        fileout = self.outputfolder + self.chromosomename + '.tab'
        try:
            SeqIO.convert(self.fastasequence, 'fasta', fileout, 'tab')
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg8)
            sys.exit(1)
        else:
            self.filelog.write(msg7)
            return fileout

    def chrcount(self):
        """
        Rename chr id of reference fasta
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
        try:
            outid = self.out + self.chromosomename + '.tab'
            self.numchr = sum(1 for line in open(self.ref2tabular()))
            self.df1 = pd.read_csv(self.ref2tabular(), header=None, sep='\t')
            self.df1[0] = self.df1.apply(lambda x: self.chromosomename + '_' + str(x.name + 1), axis=1)
            self.df1.to_csv(outid, header=None, sep='\t', index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg114)
            sys.exit(1)
        else:
            self.filelog.write(msg115)
            return outid

    def fastareference(self):
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
        try:
            with open(self.outputfolder + self.chromosomename + '.fasta', 'w') as self.ref:
                print self.ref
                for self.seq_record in SeqIO.parse(self.fastasequence, "fasta"):
                    self.ref.write('>' + str(self.chromosomename) + '\n' + str(self.seq_record.seq))
        except OSError:
            self.filelog.write(msg43)
            sys.exit(1)
        else:
            for item in self.dir:
                if item.endswith(".fai"):
                    os.remove(self.outputfolder+item)
            self.filelog.write(msg44)
            self.index = pysam.Fastafile(self.outputfolder + self.chromosomename + '.fasta')
            print self.index
            return self.outputfolder + self.chromosomename + '.fasta'

    def tab2fasta(self):
        """
        Convert tabular fasta file (2 columns) first ID second nucleotide sequence in fasta file format
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            with open(self.out + self.chromosomename + '.tab', 'r') as f:
                with open(self.out + 'reformat.fasta', 'w') as f_out:
                    for line in f:
                        line = line.strip().split('\t')
                        self.header = '>' + '_'.join([line[i] for i in self.id])
                        f_out.write(self.header + '\n')
                        f_out.write(line[self.seqix] + '\n')
                f_out.close()
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg100)
            sys.exit(1)
        else:
            self.filelog.write(msg101)
            return self.out + 'reformat.fasta'


class AnnotationFile(GenomeFile):

    def __init__(self, optparseinstance):
        GenomeFile.__init__(self, optparseinstance)
        self.NameFileAnnotation = os.path.basename(self.annotation.split('/')[-1])
        self.PathAnnotation = os.path.dirname(self.annotation)
        self.dfAnno = None
        self.aux = None
        self.fastacds = None
        self.df1 = None
        self.df2 = None

    def annotationbuild(self):
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
        try:
            self.dfAnno = pd.read_csv(self.annotation, sep="\t", header=2)
            self.aux = self.dfAnno['Location'].apply(lambda x: x.split('..'))
            self.dfAnno['Start'] = self.aux.apply(lambda x: x[0])
            self.dfAnno['End'] = self.aux.apply(lambda x: x[1])
            self.dfAnno['Chr'] = self.chromosomename
            self.dfAnno.loc[self.dfAnno['Strand'] == '+', ['Start']] = self.dfAnno['Start'].astype(int).sub(1)
            self.dfAnno[['Chr', 'Start', 'End', 'Synonym', 'COG', 'Strand', 'Gene',
                         'Product']].to_csv(self.outputfolder + self.chromosomename + '_proteome.bed',
                                            header=None, sep='\t', index=False)
        except OSError:
            self.filelog.write(msg41)
            sys.exit(0)
        else:
            self.filelog.write(msg42)
            return self.outputfolder + self.chromosomename + '_proteome.bed'








