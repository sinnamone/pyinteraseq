from output_message import *
import sys
import os
import datetime
import subprocess
from Bio import SeqIO
import gzip


class InputCheck(object):

    def __init__(self, optparseinstance):
        # import instance with all input flag
        self.inputistance = optparseinstance
        #
        self.readforward = self.inputistance.readforward
        self.readreverse = self.inputistance.readreverse
        self.outputfolder = self.inputistance.outputfolder
        self.outputid = self.inputistance.outputid
        self.primer5forward = self.inputistance.primer5forward
        self.primer3forward = self.inputistance.primer3forward
        self.primer5reverse = self.inputistance.primer5reverse
        self.primer3reverse = self.inputistance.primer3reverse
        self.fastasequence = self.inputistance.fastasequence
        self.thread = self.inputistance.thread

        #
        # self.filelog contain string of created file log
        self.minclonelength = self.inputistance.minclonelength
        self.overlapintersect = self.inputistance.overlapintersect
        self.opengap = self.inputistance.opengap
        self.mismatch = self.inputistance.mismatch
        #
        self.count = 0
        self.out = self.outputfolder + self.outputid
        self.cutadapt = '/usr/local/bin/cutadapt'
        self.trimmomatic = '/opt/Trimmomatic-0.36/trimmomatic-0.36.jar'
        self.samtools = '/usr/bin/samtools'
        self.bowtie2 = '/usr/local/bin/bowtie2'
        self.bowtie2build = '/usr/local/bin/bowtie2-build'
        self.cutadaptversion = "1.12"
        self.pick_otus = "/opt/miniconda3/envs/qiime1/bin/pick_otus.py"
        self.pick_rep_set = "/opt/miniconda3/envs/qiime1/bin/pick_rep_set.py"
        # check if input forward file is gz
        self.namefilefasta = os.path.basename(self.fastasequence.split('/')[-1])
        self.genename = self.namefilefasta.split('.')[0]
        # check file log

        #
        self.sequencingtype = self.sequencingtypecheck()
        self.first_line = None
        self.readforwardtype = self.readtypecheck(self.readforward)
        self.readreversetype = self.readtypecheck(self.readreverse)
        # assign log string
        self.inputfilelog = self.logfilecreation()


    def outputfoldercheck(self):
        if self.outputfolder is not None:
            if self.outputfolder.endswith('/') is True:
                self.outputfolder = self.outputfolder
            else:
                self.outputfolder = self.outputfolder + '/'
        else:
            sys.exit(1)
        return self.outputfolder

    def outpuedidcheck(self):
        if self.outputid is not None:
            self.outputid = self.outputid
        elif self.outputid is None:
            sys.exit(1)

    def sequencingtypecheck(self):
        if self.readreverse is None:
            self.sequencingtype = 'Single-End'
            self.readforward = self.gzipopen(self.readforward, "_forward.")
        else:
            self.sequencingtype = 'Paired-End'
            self.readforward = self.gzipopen(self.readforward, "_forward.")
            self.readreverse = self.gzipopen(self.readreverse, "_reverse.")
        return self.sequencingtype

    def readtypecheck(self, readfile):
        if readfile is not None:
            return self.fastatesting(readfile)

    def logopen(self):
        """
        Call function that open log file
        :return:
        """
        if self.inputistance.log is None:
            return open(self.outputfolder + self.outputid + "_mapping.log", "a")
        else:
            return open(self.inputfilelog, "a")

    def gzipopen(self, readfile, direction):
        """

        :param readfile:
        :param direction:
        :return:
        """
        if str(readfile).endswith(".gz"):
            outfile = readfile.split(".")
            outstring = self.outputfolder + self.outputid + direction + outfile[1]
            with gzip.open(readfile, "rt") as handle:
                with open(outstring, "w") as outfastq:
                    filecontent = handle.read()
                    outfastq.write(filecontent)
            return outstring
        else:
            return readfile

    def logfilecreation(self):
        """
        Function that open new log file.
        :return: path + name file log
        """
        if os.access(self.outputfolder, os.W_OK) is True:
            if self.inputistance.log is None:
                open(self.outputfolder + self.outputid + "_mapping.log", "a")
                return self.outputfolder + self.outputid + "_mapping.log"
            else:
                return self.inputistance.log
        else:
            sys.stdout.write(msg13)
            sys.exit(1)

    def fastatesting(self, seqinput):
        """

        :param seqinput:
        :return:
        """
        f = open(seqinput, "r")
        self.first_line = f.readline().rstrip('\n')
        if self.first_line[0] == '@':
            return "fastq"
        elif self.first_line[0] == '>':
            return "fasta"

    def checkreversereads(self, readreverse):
        """
        Check if input is paired or single end
        :return:
        """
        # check for reverse read and assign sequencing type dataset
        if readreverse is None:
            self.sequencingtype = 'Single-End'
            return self.sequencingtype
        else:
            self.sequencingtype = 'Paired-End'
            return self.sequencingtype

    def checkreads(self, varreads, message):
        """
        Chech if forward or reverse reads exist
        :param varreads:
        :param message:
        :return:
        """
        if varreads is None:
            self.logopen()
            self.filelog.write(message)
            sys.exit(1)
        else:
            return varreads

    def checkfasta(self):
        """

        :return:
        """
        # check fasta sequence
        if self.inputistance.fastasequence is None:
            self.filelog.write(msg11)
            sys.exit(1)
        else:
            return self.fastasequence

    def checkgenename(self):
        """

        :return:
        """
        self.logopen()
        # check if the name of chromosome was given by the user
        if self.namefilefasta.split('.')[0] is None:
            self.inputfilelog.write(msg87)
            sys.exit(1)
        else:
            return self.genename

    def fastqcount(self, fastq, rtype):
        """
        Function to count the number of sequence
        :param fastq:
        :param rtype:
        :return:
        """
        self.count = 0
        for record in SeqIO.parse(fastq, rtype):
            self.count = self.count + 1
        return str(self.count)

    def inputinformationappen(self):
        """
        Log compilation. Call all the previous functions
        :return:
        """
        self.logopen().write(datetime.datetime.now().ctime() + '\n')
        self.logopen().write(msg14)
        self.logopen().write(msg15 + self.checkreversereads(readreverse=self.readreverse))

        if self.sequencingtype == 'Paired-End':
            # check dataset
            self.logopen().write(msg17 + self.checkreads(varreads=self.readforward,
                                                           message=msg2))
            self.logopen().write(msg18 + self.checkreads(varreads=self.readreverse,
                                                           message=msg3))
            self.logopen().write(msg21 + self.primer5forward)
            self.logopen().write(msg22 + self.primer3forward)
            self.logopen().write(msg23 + self.primer5reverse)
            self.logopen().write(msg24 + self.primer3reverse)
            # self.openlogfile.write(msg47 + self.fastqcount(fastq=self.readforward,
            #                                                rtype=self.readforwardtype))
            # self.openlogfile.write(msg48 + self.fastqcount(fastq=self.readreverse,
            #                                                rtype=self.readreversetype)) #TODO
        else:
            self.readforward = self.gzipopen(self.readforward, "_forward.")
            self.logopen().write(msg28 + self.checkreads(varreads=self.readforward,
                                                         message=msg1))
            self.logopen().write(msg29 + self.readforwardtype)
            self.logopen().write(msg30 + self.primer5forward)
            self.logopen().write(msg31 + self.primer3forward)
            self.logopen().write(msg49 + self.fastqcount(fastq=self.readforward,
                                                         rtype=self.readforwardtype))
        self.logopen().write(msg25 + self.outputid)
        self.logopen().write(msg26 + self.checkfasta())
        self.logopen().write(msg27 + self.checkgenename())
        self.logopen().write(msg33 + subprocess.check_output(['cutadapt', '--version']))
        self.logopen().write(msg0)
        self.logopen().close()
        return True
