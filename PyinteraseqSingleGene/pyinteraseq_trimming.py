from output_message import *
from pyinteraseq_inputcheck import InputCheck
import sys
import subprocess


class TrimmingSingle(InputCheck):

    def __init__(self, optparseinstance):
        InputCheck.__init__(self, optparseinstance)
        self.primer5reverse = InputCheck(optparseinstance)

    def trimming5single(self):
        """
        Function for trimming 5' single-end reads
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_trimming.log", "a")
        try:
            subprocess.check_call(['cutadapt', '-g', str(self.primer5forward), '--discard-untrimmed',
                                   '-e', '0.03', '--trim-n', '-m', str(self.cloneslength), '--quiet',
                                   '-o', self.out + '_read1.' + self.readforwardtype,
                                   str(self.readforward)], stderr=self.filelog)
        except subprocess.CalledProcessError:
            self.filelog.write(msg34)
            sys.exit(1)
        else:
            self.filelog.write(msg50 + self.fastqcount(self.out + '_read1.' + self.readforwardtype,
                                                       self.readforwardtype))
            self.filelog.write(msg35)
            return self.out + '_read1.' + self.readforwardtype

    def trimming3single(self):
        """
        Function for trimming 3' single-end reads
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_trimming.log", "a")
        try:
            subprocess.check_call(
                ['cutadapt', '-a', str(self.primer3forward), '--trim-n',
                 '-e', '0.03', '-m', str(self.cloneslength), '--quiet',
                 '-o', self.out + '_read1_1.' + self.readforwardtype,
                 self.out + '_read1.' + self.readreversetype],
                stderr=self.filelog)
        except subprocess.CalledProcessError:
            self.filelog.write(msg36)
            sys.exit(1)
        else:
            self.filelog.write(msg51 + self.fastqcount(self.out + '_read1_1.' + self.readforwardtype,
                                                       self.readforwardtype))
            self.filelog.write(msg37)
            return self.out + '_read1_1.' + self.readforwardtype


class TrimmingPaired(InputCheck):

    def __init__(self, optparseinstance):
        InputCheck.__init__(self, optparseinstance)
        self.readtrimmedfive = []
        self.readtrimmedthree = []

    def trimming5paired(self):
        """
        Trimming 5' paired-end dataset
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_trimming.log", "a")
        try:
            subprocess.check_call(['cutadapt',
                                   '-g', self.primer5forward,
                                   '-G', self.primer5reverse,
                                   '-e', '0.03',
                                   '--trim-n',
                                   '-m', self.cloneslength,
                                   '--quiet',
                                   '--discard-untrimmed',
                                   '-o', self.out + '_read1.' + self.readforwardtype,
                                   '-p', self.out + '_read2.' + self.readreversetype,
                                   self.readforward,
                                   self.readreverse],
                                  stderr=self.filelog)
        except subprocess.CalledProcessError:
            self.filelog.write(msg38a)
            sys.exit(1)
        else:
            self.readtrimmedfive.append(self.out + '_read1.' + self.readforwardtype)
            self.readtrimmedfive.append(self.out + '_read2.' + self.readreversetype)
            self.filelog.write(msg52 + self.fastqcount(self.out + '_read1.' + self.readforwardtype,
                                                       self.readforwardtype))
            self.filelog.write(msg53 + self.fastqcount(self.out + '_read2.' + self.readreversetype,
                                                       self.readreversetype))
            self.filelog.write(msg38b)
            return self.readtrimmedfive

    def trimming3paired(self):
        """
        Trimming 3' paired-end
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_trimming.log", "a")
        try:
            subprocess.check_call(
                ['cutadapt',
                 '-b', self.primer3forward,
                 '-B', self.primer3reverse,
                 '-e', '0.03',
                 '--trim-n', '-m', self.cloneslength,
                 '--quiet',
                 '-o', self.out + '_read1_1.' + self.readforwardtype,
                 '-p', self.out + '_read2_2.' + self.readreversetype,
                 self.out + '_read1.' + self.readforwardtype,
                 self.out + '_read2.' + self.readreversetype],
                stderr=self.filelog)
        except subprocess.CalledProcessError:
            sys.stdout.write(msg39)
            sys.exit(1)
        else:
            self.readtrimmedthree.append(self.out + '_read1_1.' + self.readforwardtype)
            self.readtrimmedthree.append(self.out + '_read2_2.' + self.readforwardtype)
            self.filelog.write(msg54 + self.fastqcount(self.out + '_read1_1.' + self.readforwardtype,
                                                       self.readforwardtype))
            self.filelog.write(msg55 + self.fastqcount(self.out + '_read2_2.' + self.readreversetype,
                                                       self.readreversetype))
            self.filelog.write(msg40)
        return self.readtrimmedthree
