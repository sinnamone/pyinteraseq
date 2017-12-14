from output_message import *
from pyinteraseq_inputcheck import InputCheck
import sys
import subprocess


class Trimming(InputCheck):

    def __init__(self, optparseinstance):
        InputCheck.__init__(self, optparseinstance)
        self.readtrimmedfive = []
        self.readtrimmedthree = []
        self.filelog = None

    def trimming5single(self, readfile, primer5, direction, readtype):
        """
        Function for trimming 5' single-end reads
        :return:
        """
        try:
            subprocess.check_call(['cutadapt', '-b', primer5,
                                   '--discard-untrimmed',
                                   '--trim-n',
                                   '-m', self.minclonelength,
                                   '--quiet',
                                   '-o', self.out + direction + readtype,
                                   readfile])
        except subprocess.CalledProcessError:
            self.logopen().write(msg34)
            sys.exit(1)
        else:
            # self.logopen().write(msg50 + self.fastqcount(self.out + direction + readtype, readtype))
            self.logopen().write(msg35)
            return self.out + direction + readtype

    def trimming3single(self, readfile, primer3, direction, readtype):
        """
        Function for trimming 3' single-end reads
        :return:
        """
        self.filelog = self.logopen()
        try:
            subprocess.check_call(['cutadapt', '-a', primer3,
                                   '--trim-n',
                                   '-m', self.minclonelength,
                                   '--quiet',
                                   '-o', self.out + direction + readtype,
                                   readfile])
        except subprocess.CalledProcessError:
            self.filelog.write(msg36)
            sys.exit(1)
        else:
            # self.filelog.write(msg51 + self.fastqcount(self.out + direction + readtype,
            #                                            readtype))
            self.filelog.write(msg37)
            return self.out + direction + readtype

    def trimming5paired(self):
        """
        Trimming 5' paired-end dataset
        :return:
        """
        try:
            subprocess.check_call(['cutadapt',
                                   '-g', self.primer5forward,
                                   '-G', self.primer5reverse,
                                   '--trim-n',
                                   '-m', self.minclonelength,
                                   '--quiet',
                                   '--discard-untrimmed',
                                   '-o', self.out + '_read1.' + self.readforwardtype,
                                   '-p', self.out + '_read2.' + self.readreversetype,
                                   self.readforward,
                                   self.readreverse])
        except subprocess.CalledProcessError:
            self.logopen().write(msg38a)
            sys.exit(1)
        else:
            self.readtrimmedfive.append(self.out + '_read1.' + self.readforwardtype)
            self.readtrimmedfive.append(self.out + '_read2.' + self.readreversetype)
            self.logopen().write(msg52 + self.fastqcount(self.out + '_read1.' + self.readforwardtype,
                                                         self.readforwardtype))
            self.logopen().write(msg53 + self.fastqcount(self.out + '_read2.' + self.readreversetype,
                                                         self.readreversetype))
            self.logopen().write(msg38b)
            return self.readtrimmedfive

    def trimming3paired(self):
        """
        Trimming 3' paired-end
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
        try:
            subprocess.check_call(
                ['cutadapt',
                 '-b', self.primer3forward,
                 '-B', self.primer3reverse,
                 '-e', '0.03',
                 '--trim-n', '-m', self.minclonelength,
                 '--quiet',
                 '-o', self.out + '_read1_1.' + self.readforwardtype,
                 '-p', self.out + '_read2_2.' + self.readreversetype,
                 self.out + '_read1.' + self.readforwardtype,
                 self.out + '_read2.' + self.readreversetype],
                stderr=self.filelog)
        except subprocess.CalledProcessError:
            self.filelog.write(msg39)
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
