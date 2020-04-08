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

    def concatenateforrev(self, readlist):
        """
        Merge forward and reverse fasta file
        :param readlist: list with files to append
        :return: path + idanalysis + _con.fasta
        """
        self.filelog = self.logopen()
        try:
            with open(self.out + '_con.' + self.readforwardtype, 'w') as outfile:
                for fname in readlist:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
            outfile.close()
        except StandardError:
            self.filelog.write(msg88)
            sys.exit(1)
        else:
            self.filelog.write(msg56 + self.fastqcount(self.out + '_con.' + self.readforwardtype,
                                                       self.readforwardtype))
            return self.out + '_con.' + self.readforwardtype
