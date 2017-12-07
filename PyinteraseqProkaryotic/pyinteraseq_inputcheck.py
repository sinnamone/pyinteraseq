from output_message import *
import sys
import os
import datetime
import subprocess
from Bio import SeqIO
import traceback
import gzip


class InputCheck(object):

    def __init__(self, optparseinstance):
        # import instance with all input flag
        self.inputistance = optparseinstance
        # put all under this
        self.count = 0
        self.inputfilelog = self.inputistance.log
        self.cutadapt = ''
        self.cutadaptversion = ""
        self.readforward = self.inputistance.readforward
        self.readreverse = self.inputistance.readreverse
        if self.inputistance.outputid is not None:
            self.outputid = self.inputistance.outputid
        if self.inputistance.outputfolder is not None:
            if self.inputistance.outputfolder.endswith('/') is True:
                self.outputfolder = self.inputistance.outputfolder
                self.out = self.outputfolder + self.outputid
            else:
                self.outputfolder = self.inputistance.outputfolder + '/'
                self.out = self.outputfolder + self.outputid
        else:
            sys.exit(1)
        if self.readreverse is None:
            self.sequencingtype = 'Single-End'
            self.readforward = self.gzipopen(self.readforward, "_forward.")
        else:
            self.sequencingtype = 'Paired-End'
            self.readforward = self.gzipopen(self.readforward, "_forward.")
            self.readreverse = self.gzipopen(self.readreverse, "_reverse.")

        self.primer5forward = self.inputistance.primer5forward
        self.primer3forward = self.inputistance.primer3forward
        self.primer5reverse = self.inputistance.primer5reverse
        self.primer3reverse = self.inputistance.primer3reverse
        self.fastasequence = self.inputistance.fastasequence
        self.pick_otus = "/opt/miniconda3/envs/qiime1/bin/pick_otus.py"
        self.pick_rep_set = "/opt/miniconda3/envs/qiime1/bin/pick_rep_set.py"
        # check if input forward file is gz
        self.thread = self.inputistance.thread
        self.cloneslength = self.inputistance.minclonelength
        self.namefilefasta = os.path.basename(self.fastasequence.split('/')[-1])
        self.genename = self.namefilefasta.split('.')[0]
        self.opengap = self.inputistance.opengap
        self.mismatch = self.inputistance.mismatch
        # check file log
        self.inputfilelog = self.logfilecreation()
        self.filelog = None
        #
        self.readforwardtype = self.fastatesting(self.readforward)
        if self.readreverse is not None:
            self.readreversetype = self.fastatesting(self.readreverse)

    def logopen(self):
        if self.inputfilelog is None:
            return open(self.outputfolder + self.outputid + "_mapping.log", "a")
        else:
            return open(str(self.inputfilelog), 'a')

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
            if self.inputfilelog is None:
                self.inputfilelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
                return self.outputfolder + self.outputid + "_mapping.log"
            else:
                return self.inputfilelog
        else:
            sys.stdout.write(msg13)
            sys.exit(1)


    def fastatesting(self, seqinput):
        """

        :param seqinput:
        :return:
        """
        f = open(seqinput, "r")
        first_line = f.readline().rstrip('\n')
        if first_line[0] == '@':
            return "fastq"
        elif first_line[0] == '>':
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
            self.readreverse = self.inputistance.readreverse
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
        self.filelog = open(self.inputfilelog, 'a')
        # check if the name of chromosome was given by the user
        if self.namefilefasta.split('.')[0] is None:
            self.filelog.write(msg87)
            sys.exit(1)
        else:
            return self.genename

    def cutadaptcheck(self):
        """
        Function that verify the installation of cutadapt
        :return: folder where cutadapt is located.
        """
        try:
            self.cutadapt = subprocess.check_output(['which', 'cutadapt']).split('\n')[0]
        except subprocess.CalledProcessError:
            self.filelog.write(msg32)
            sys.exit(1)
        else:
            return self.cutadapt

    def pickotuscheck(self):
        """
        Function that verify the installation of pick_otus
        :return:
        """
        try:
            self.pick_otus = subprocess.check_output(['which', 'pick_otus.py']).split('\n')[0]
        except subprocess.CalledProcessError:
            self.filelog.write(msg67)
            sys.exit(1)
        else:
            return self.pick_otus

    def pickrepseqcheck(self):
        """
        Function that verify the installation of pick_rep_set
        :return:
        """
        try:
            self.pick_rep_set = subprocess.check_output(['which', 'pick_rep_set.py']).split('\n')[0]
        except subprocess.CalledProcessError:
            self.filelog.write(msg68)
            sys.exit(1)
        else:
            return self.pick_rep_set

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
        self.filelog = self.logopen()
        self.filelog.write(datetime.datetime.now().ctime() + '\n')
        self.filelog.write(msg14)
        self.filelog.write(msg15 + self.checkreversereads(readreverse=self.readreverse))

        if self.sequencingtype == 'Paired-End':
            # check dataset
            self.readforward = self.gzipopen(self.readforward, "_forward.")
            self.readreverse = self.gzipopen(self.readreverse, "_reverse.")
            self.filelog.write(msg17 + self.checkreads(varreads=self.readforward,
                                                       message=msg2))
            self.filelog.write(msg18 + self.checkreads(varreads=self.readreverse,
                                                       message=msg3))
            # check primers
            self.filelog.write(msg21 + self.primer5forward)
            self.filelog.write(msg22 + self.primer3forward)
            self.filelog.write(msg23 + self.primer5reverse)
            self.filelog.write(msg24 + self.primer3reverse)
            self.filelog.write(msg47 + self.fastqcount(fastq=self.readforward,
                                                       rtype=self.readforwardtype))
            self.filelog.write(msg48 + self.fastqcount(fastq=self.readreverse,
                                                       rtype=self.readreversetype))
        else:
            self.readforward = self.gzipopen(self.readforward, "_forward.")
            self.filelog.write(msg28 + self.checkreads(varreads=self.readforward,
                                                       message=msg1))
            self.filelog.write(msg29 + self.readforwardtype)
            self.filelog.write(msg30 + self.primer5forward)
            self.filelog.write(msg31 + self.primer3forward)
            self.filelog.write(msg49 + self.fastqcount(fastq=self.readforward,
                                                       rtype=self.readforwardtype))
        self.filelog.write(msg25 + self.outputid)
        self.filelog.write(msg26 + self.checkfasta())
        self.filelog.write(msg27 + self.checkgenename())
        self.filelog.write(msg33 + subprocess.check_output(['cutadapt', '--version']))
        self.filelog.write(msg0)
        self.filelog.close()
        return True


class InputCheckDomainDefinition(object):

        def __init__(self, optparseinstance):
            # import instance with all input flag
            self.inputistance = optparseinstance
            # put all under this
            self.thread = self.inputistance.thread
            self.count = 0
            self.inputfilelog = self.inputistance.log
            self.cloneslength = self.inputistance.minclonelength
            self.outputfolder = self.checkoutputpath()
            self.outputid = self.checkoutputid()
            self.fastasequence = self.checkfastasequence()
            self.namefilefasta = os.path.basename(self.fastasequence.split('/')[-1])
            self.genename = self.namefilefasta.split('.')[0]
            self.out = self.outputfolder + self.outputid
            self.pick_otus = "/opt/miniconda3/envs/qiime1/bin/pick_otus.py"
            self.pick_rep_set = "/opt/miniconda3/envs/qiime1/bin/pick_rep_set.py"
            # check file log
            self.inputfilelog = self.logfilecreation()
            self.filelog = None

            # check annotation
            if self.inputistance.annotation is not None:
                self.annotation = self.inputistance.annotation

        def logfilecreation(self):
            """
            Function that open new log file.
            :return: path + name file log
            """
            if os.access(self.outputfolder, os.W_OK) is True:
                if self.inputfilelog is None:
                    self.inputfilelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
                    return self.outputfolder + self.outputid + "_mapping.log"
                else:
                    return self.inputfilelog
            else:
                sys.stdout.write(msg13)
                sys.exit(1)

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

        def checkoutputpath(self):
            """

            :return:
            """
            if self.inputistance.outputfolder is not None:
                if self.inputistance.outputfolder.endswith('/') is True:
                    self.outputfolder = self.inputistance.outputfolder
                else:
                    self.outputfolder = self.inputistance.outputfolder + '/'
                return self.outputfolder

        def checkoutputid(self):
            """

            :return:
            """
            if self.inputistance.outputid is None:
                self.filelog.write(msg10)
                sys.exit(1)
            else:
                return self.inputistance.outputid

        def logopen(self):
            if self.inputfilelog is None:
                return open(self.outputfolder + self.outputid + "_domaind_definition.log", "a")
            else:
                return open(str(self.inputfilelog), 'a')

        def checkfastasequence(self):
            if self.inputistance.fastasequence is None:
                self.filelog.write(msg11)
                sys.exit(1)
            else:
                return self.inputistance.fastasequence


