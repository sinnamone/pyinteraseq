import subprocess
from output_message import *
import os
import sys
from pyinteraseq_trimming import Trimming
from multiprocessing import Pool
import traceback
import urllib
import collections
import pysam
from Bio import SeqIO
import numpy as np


class BlastNlucleotide(Trimming):

    def __init__(self, optparseinstance):
        Trimming.__init__(self, optparseinstance)
        self.header = []
        self.pool = Pool(processes=int(self.thread))
        self.dbname = self.outputfolder + os.path.basename(self.fastasequence.split('/')[-1]).split('.')[0]
        self.filelog = None
        self.parts = ''
        self.samInfoFields = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL",
                              "ASCORE", "ABASES", "MISMATCHES", "GAPOPENS", "GAPEXTENSIONS", "EDITDISTANCE",
                              "MISMATCHEDREFBASES", "PAIR"]
        self.SAMRecord = collections.namedtuple("SAMRecord", self.samInfoFields)
        self.samInfoFallistoFields = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ","QUAL","NH"]
        self.SAMRecordKallisto = collections.namedtuple("SAMRecord", self.samInfoFallistoFields)

    def indexing_bowtie(self):
        """

        :return:
        """
        self.filelog = self.logopen()
        fnull = open(os.devnull, 'w')
        try:
            if os.path.exists(self.dbname + '.1.bt2') and os.path.exists(self.dbname + '.rev.1.bt2'):
                self.filelog.write(msg123)
            else:
                subprocess.check_call(['/usr/local/bin/bowtie2-build', "--threads", self.thread,
                                       '-q', '-f', self.fastasequence,
                                       self.dbname],
                                      stderr=fnull, stdout=fnull)
        except subprocess.CalledProcessError:
            self.filelog.write(msg43)
            sys.exit(1)
        except traceback:
            self.filelog.write(msg43)
            sys.exit(1)
        else:
            self.filelog.write(msg44)
            return self.dbname

    def mapping_bowtie(self, multifasta, database):
        """
        Function that call bowtie2 for mapping.
        :param multifasta:
        :param database:
        :return:
        """
        self.filelog = self.logopen()
        fnull = open(os.devnull, 'w')
        try:
            if self.readforwardtype == 'fastq':
                subprocess.check_call(
                    ["/usr/local/bin/bowtie2", '--fast-local', '-p', '20', '-x',
                     database, '-q', multifasta, ' -S', self.out + '.sam'],
                    stderr=fnull, stdout=fnull)
            elif self.readforwardtype == 'fasta':
                subprocess.check_call(
                    ["/usr/local/bin/bowtie2", '--fast-local', '-p', self.thread, '-x',
                     database, '-f', multifasta, ' -S', self.out + '.sam'],
                    stderr=fnull, stdout=fnull)
        except subprocess.CalledProcessError:
            self.filelog.write(msg60)
            sys.exit(1)
        else:
            self.filelog.write(msg61)
            return self.out + '.sam'

    def getheadersam(self, samfile):
        """
        Get header Sam file
        :param samfile:
        :return:
        """
        self.filelog = self.logopen()
        try:
            with open(samfile) as infile:
                for line in infile:
                    if line.startswith("@"):
                        self.header.append(line)
        except traceback:
            self.filelog.write(msg58)
            sys.exit(1)
        else:
            self.filelog.write(msg59)
            return self.header

    def parsesamfiles(self, samfile):
        """

        :param samfile:
        :return:
        """
        self.filelog = self.logopen()
        try:
            with open(samfile) as infile:
                for line in infile:
                    if line.startswith("@"):
                        continue
                    if line.split("\t")[5] is '*':
                        continue
                    if line.split("\t")[12].startswith('XS:i:'):
                        continue
                    self.parts = line.split("\t")
                    assert len(self.parts) == len(self.samInfoFields)
                    # Normalize data
                    normalizedinfo = {
                        "QNAME": None if self.parts[0] == "." else urllib.unquote(self.parts[0]),
                        "FLAG": None if self.parts[1] == "." else urllib.unquote(self.parts[1]),
                        "RNAME": None if self.parts[2] == "." else urllib.unquote(self.parts[2]),
                        "POS": None if self.parts[3] == "." else urllib.unquote(self.parts[3]),
                        "MAPQ": None if self.parts[4] == "." else urllib.unquote(self.parts[4]),
                        "CIGAR": None if self.parts[5] == "." else urllib.unquote(self.parts[5]),
                        "RNEXT": None if self.parts[6] == "." else urllib.unquote(self.parts[6]),
                        "PNEXT": None if self.parts[7] == "." else urllib.unquote(self.parts[7]),
                        "TLEN": None if self.parts[8] == "." else urllib.unquote(self.parts[8]),
                        "SEQ": None if self.parts[9] == "." else urllib.unquote(self.parts[9]),
                        "QUAL": None if self.parts[10] == "." else urllib.unquote(self.parts[10]),
                        "ASCORE": None if self.parts[11] == "." else urllib.unquote(self.parts[11]),
                        "ABASES": None if self.parts[12] == "." else urllib.unquote(self.parts[12]),
                        "MISMATCHES": None if self.parts[13] == "." else urllib.unquote(self.parts[13]),
                        "GAPOPENS": None if self.parts[14] == "." else urllib.unquote(self.parts[14]),
                        "GAPEXTENSIONS": None if self.parts[15] == "." else urllib.unquote(self.parts[11]),
                        "EDITDISTANCE": None if self.parts[16] == "." else urllib.unquote(self.parts[12]),
                        "MISMATCHEDREFBASES": None if self.parts[17] == "." else urllib.unquote(self.parts[13]),
                        "PAIR": None if self.parts[18] == "." else urllib.unquote(self.parts[14])
                    }
                    yield self.SAMRecord(**normalizedinfo)
        except traceback:
            self.filelog.write(msg63)
            sys.exit(1)
        else:
            self.filelog.write(msg64)

    def parsesamkallistofiles(self, samfile):
        """

        :param samfile:
        :return:
        """
        self.filelog = self.logopen()
        try:
            with open(samfile) as infile:
                for line in infile:
                    if line.startswith("@"):
                        continue
                    if line.split("\t")[5] is '*':
                        continue
                    self.parts = line.split("\t")
                    assert len(self.parts) == len(self.samInfoFallistoFields)
                    # Normalize data
                    normalizedinfo = {
                        "QNAME": None if self.parts[0] == "." else urllib.unquote(self.parts[0]),
                        "FLAG": None if self.parts[1] == "." else urllib.unquote(self.parts[1]),
                        "RNAME": None if self.parts[2] == "." else urllib.unquote(self.parts[2]),
                        "POS": None if self.parts[3] == "." else urllib.unquote(self.parts[3]),
                        "MAPQ": None if self.parts[4] == "." else urllib.unquote(self.parts[4]),
                        "CIGAR": None if self.parts[5] == "." else urllib.unquote(self.parts[5]),
                        "RNEXT": None if self.parts[6] == "." else urllib.unquote(self.parts[6]),
                        "PNEXT": None if self.parts[7] == "." else urllib.unquote(self.parts[7]),
                        "TLEN": None if self.parts[8] == "." else urllib.unquote(self.parts[8]),
                        "SEQ": None if self.parts[9] == "." else urllib.unquote(self.parts[9]),
                        "QUAL": None if self.parts[10] == "." else urllib.unquote(self.parts[10]),
                        "NH": None if self.parts[10] == "." else urllib.unquote(self.parts[11])
                    }
                    yield self.SAMRecordKallisto(**normalizedinfo)
        except traceback:
            self.filelog.write(msg63)
            sys.exit(1)
        else:
            self.filelog.write(msg64)

    def filteringmismatches(self):
        self.filelog = self.logopen()
        try:
            target = open(self.out + '_filtered.sam', 'w')
            for i in range(len(self.header)):
                target.write(self.header[i].strip('\n') + '\n')
            for record in self.parsesamfiles(samfile=self.out + '.sam'):
                if (float(record.MISMATCHES.split(':')[-1]) * 100) / len(record.SEQ) <= self.mismatch:
                    target.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tAS:i:%s\n" % (record.QNAME,
                                                                                            record.FLAG,
                                                                                            record.RNAME,
                                                                                            record.POS,
                                                                                            record.MAPQ,
                                                                                            record.CIGAR,
                                                                                            record.RNEXT,
                                                                                            record.PNEXT,
                                                                                            record.TLEN,
                                                                                            record.SEQ,
                                                                                            record.QUAL,
                                                                                            (int(
                                                                                                record.MISMATCHES.split(
                                                                                                    ':')[
                                                                                                    -1]) * 100) / len(
                                                                                                record.SEQ)))

            target.close()
        except traceback:
            self.filelog.write(msg80)
            sys.exit(1)
        else:
            self.filelog.write(msg81)
            return self.out + '_filtered.sam'

    def filteringmismatcheskallisto(self):
        self.filelog = self.logopen()
        try:
            target = open(self.out + '_filtered.sam', 'w')
            for i in range(len(self.header)):
                target.write(self.header[i].strip('\n') + '\n')
            for record in self.parsesamkallistofiles(samfile=self.out + '.sam'):
                target.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (record.QNAME,
                                                                                   record.FLAG,
                                                                                   record.RNAME,
                                                                                   record.POS,
                                                                                   record.MAPQ,
                                                                                   record.CIGAR,
                                                                                   record.RNEXT,
                                                                                   record.PNEXT,
                                                                                   record.TLEN,
                                                                                   record.SEQ,
                                                                                   record.QUAL,
                                                                                   record.NH))
            target.close()
        except traceback:
            self.filelog.write(msg80)
            sys.exit(1)
        else:
            self.filelog.write(msg81)
            return self.out + '_filtered.sam'

    def conversionsam2bam(self, samfile):
        self.filelog = self.logopen()
        try:
            pysam.view("-hbS", samfile, "-F", "4", "-@", self.thread, "-o", self.out + ".bam", catch_stdout=False)
        except traceback:
            self.filelog.write(msg82)
            sys.exit(1)
        else:
            self.filelog.write(msg83)
            return self.out + ".bam"

    def sortbam(self, bamfile):
        """
        Sort and indexing of BAM file
        :param bamfile:
        :return:
        """
        self.filelog = self.logopen()
        try:
            pysam.sort("-@", self.thread, bamfile, "-o", self.out + "_mapping.bam")
            pysam.index(self.out + "_mapping.bam")
        except subprocess.CalledProcessError:
            self.filelog.write(msg84)
            sys.exit(1)
        else:
            self.filelog.write(msg85)
            return self.out + "_mapping.bam"

    def kallistoindexvalues(self, fastaconcatenate):
        try:
            a = []
            for seq_record in SeqIO.parse(fastaconcatenate, "fastq"):
                a.append((len(seq_record)))
        except traceback:
            self.filelog.write(msg131)
            sys.exit(1)
        else:
            self.filelog.write(msg132)
            return np.mean(a), np.std(a)

    def indexkallisto(self):
        self.filelog = self.logopen()
        fnull = open(os.devnull, 'w')
        try:
            if os.path.exists(self.dbname):
                self.filelog.write(msg123)
            else:
                subprocess.check_call(["/opt/miniconda3/bin/kallisto", "index", "-i", self.dbname, self.fastasequence],
                                      stderr=fnull, stdout=fnull)
        except subprocess.CalledProcessError:
            self.filelog.write(msg43)
            sys.exit(1)
        except traceback:
            self.filelog.write(msg43)
            sys.exit(1)
        else:
            self.filelog.write(msg44)
            return self.dbname

    def mappingkallisto(self, indexkallisto, avefragmlen, stdfragmlen, reads):
        self.filelog = self.logopen()
        fnull = open(os.devnull, 'w')
        try:
            with open(self.out + '.sam', 'w') as sam:
                subprocess.check_call(
                    ["/opt/miniconda3/bin/kallisto", 'quant', '-i', indexkallisto, '-o', self.outputfolder,
                     '--pseudobam', '--single', '-l', str(avefragmlen), '-s', str(stdfragmlen), reads],
                    stderr=fnull, stdout=sam)
        except subprocess.CalledProcessError:
            self.filelog.write(msg60)
            sys.exit(1)
        else:
            self.filelog.write(msg61)
            return self.out + '.sam'

    def cleantempfile(self):
        """
        Remove temporany files.
        :return:
        """
        templistfilesingle = [".sam", "_filtered.sam", ".bam", "_forward.fastq", "_forward5trimmed.fastq",
                              "_forward3trimmed.fastq", "_reverse3trimmed.fasta", "_forward3trimmed.fasta",
                              "_reverse5trimmed.fasta", "_forward5trimmed.fasta", "_con.fasta"]
        templistfilepaired = ["_forward.fastq", "_reverse.fastq", "_forward5trimmed.fastq", "_reverse5trimmed.fastq",
                              "_forward3trimmed.fastq", "_reverse3trimmed.fastq", "_con.fastq", ".sam", "_filtered.sam",
                              ".bam"]
        if self.sequencingtype in "Single-End":
            if self.readforwardtype in "fastq":
                for item in templistfilesingle:
                    if os.path.isfile(self.out + item):
                        os.remove(self.out + item)
            elif self.readforwardtype in "fasta":
                for item in templistfilesingle:
                    if os.path.isfile(self.out + item):
                        os.remove(self.out + item)
        elif self.sequencingtype in "Paired-End":
            if self.readforwardtype in "fastq":
                for item in templistfilepaired:
                    if os.path.isfile(self.out + item):
                        os.remove(self.out + item)
            elif self.readforwardtype in "fasta":
                for item in templistfilepaired:
                    if os.path.isfile(self.out + item):
                        os.remove(self.out + item)

