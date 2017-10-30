import os
from Bio import SeqIO
from Bio.Alphabet import ProteinAlphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import pandas as pd
import sys
from output_message import *
import pysam
import pybedtools
import warnings
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
            # save the name of chromosome. Format string is composed by >gi|XXXXXX|ref|NC_XXXXX.1|description
            # self.chromosomename = self.seq_record.id.split("|")[-2]
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
        self.fastacds = None
        self.df1 = None
        self.df2 = None

    def annotationbuild(self):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
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


# class ProteomeParsing(AnnotationFile):
#
#     def __init__(self, optparseinstance):
#         AnnotationFile.__init__(self, optparseinstance)
#         warnings.filterwarnings("ignore")
#
#     def translatednaframescds(self, seq, inputfile):
#         input_handle = open(inputfile, "r")
#         output_handle = open(self.outputfolder + self.chromosomename + '_allframescdsfasta.tab', "a+")
#         allpossibilities = []
#         for frame in range(3):
#             trans = str(seq.seq[frame:].translate(11))
#             allpossibilities.append(trans)
#         for frame in range(3):
#             trans = str(seq.seq.reverse_complement()[frame:].translate(11))
#             allpossibilities.append(trans)
#         i = 0
#         for currentFrame in allpossibilities:
#             i = i + 1
#             currentprotein = Seq(currentFrame, alphabet=ProteinAlphabet)
#
#             currentproteinrecord = SeqRecord(currentprotein, seq.name)
#             currentproteinrecord.id = currentproteinrecord.id + "." + str(i)
#             currentproteinrecord.description = seq.description + "; frame " + str(i)
#             SeqIO.write(currentproteinrecord, output_handle, "tab")
#         input_handle.close()
#         output_handle.close()
#         return self.outputfolder + self.chromosomename + '_allframescdsfasta.tab'
#
#     def getcdsfasta(self):
#         annotation = pybedtools.BedTool(self.annotationbuild())
#         fastasequence = pybedtools.BedTool(self.fastareference())
#         self.fastacds = annotation.sequence(fi=fastasequence, s=True).save_seqs(
#             self.outputfolder + self.chromosomename + '_fastacds.fasta')
#         for seq_record in SeqIO.parse(self.outputfolder + self.chromosomename + '_fastacds.fasta', "fasta",
#                                       alphabet=IUPAC.ambiguous_dna):
#             self.translatednaframescds(seq_record, self.outputfolder + self.chromosomename + '_fastacds.fasta')
#         return self.outputfolder + self.chromosomename + '_allframescdsfasta.tab'
#
#     def filtercdsframe(self):
#         self.df1 = pd.read_csv(self.getcdsfasta(), sep="\t", header=None, names=['cdsid', 'aaseq'])
#         # self.df2 = self.df1[self.df1['aaseq'].str.contains('\*') == False]
#         self.df2 = self.df1[self.df1['aaseq'].str.startswith("*") == True]
#         self.df3 = self.df1[self.df1['aaseq'].str.endswith("*")) == True]
#         self.df2.to_csv(self.outputfolder + self.chromosomename + '_allcdsframesfilteredstart.tab', sep="\t",
#                         header=None, index=False)
#         self.df3.to_csv(self.outputfolder + self.chromosomename + '_allcdsframesfilteredend.tab', sep="\t",
#                         header=None, index=False)
#         return self.outputfolder + self.chromosomename + '_allcdsframesfiltered.tab'



