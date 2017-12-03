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
        # put all under this
        # dataset input
        self.readforward = self.inputistance.readforward
        self.readreverse = self.inputistance.readreverse
        # Thread
        self.thread = self.inputistance.thread
        # Output folder
        self.outputfolder = self.checkoutputpath()
        # Output Id
        self.outputid = self.checkoutputid()
        #
        self.out = self.outputfolder + self.outputid
        # pick otus
        self.pick_otus = None
        self.pick_rep_set = None
        self.cloneslength = self.inputistance.minclonelength
        self.fastasequence = self.checkfastasequence()
        self.namefilefasta = os.path.basename(self.fastasequence.split('/')[-1])
        self.genename = self.namefilefasta.split('.')[0]
        self.opengap = self.inputistance.opengap
        self.mismatch = self.inputistance.mismatch
        self.count = 0
        self.cutadapt = ''
        self.cutadaptversion = ""
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
        self.sequencingtype = self.checksequencingtype()
        # check if input forward file is gz
        self.readforward = self.gzipopen(self.readforward, "_forward.")
        self.readreverse = self.gzipopen(self.readreverse, "_reverse.")
        # check sequencing type
        if self.sequencingtype in "Single-End":
            self.readforwardtype = self.fastatesting(self.readforward)
        elif self.sequencingtype in "Paired-End":
            self.readforwardtype = self.fastatesting(self.readforward)
            self.readreversetype = self.fastatesting(self.readreverse)
        self.primer5forward = self.inputistance.primer5forward
        self.primer3forward = self.inputistance.primer3forward
        self.primer5reverse = self.inputistance.primer5reverse
        self.primer3reverse = self.inputistance.primer3reverse

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

    def openlog(self):
        return open(self.outputfolder + self.outputid + "_mapping.log", "a")

    def checkoutputid(self):
        if self.inputistance.outputid is None:
            self.filelog.write(msg10)
            sys.exit(1)
        else:
            return self.inputistance.outputid

    def checkfastasequence(self):
        if self.inputistance.fastasequence is None:
            self.filelog.write(msg11)
            sys.exit(1)
        else:
            return self.inputistance.fastasequence

    def checksequencingtype(self):
        if self.readreverse is None:
            self.sequencingtype = 'Single-End'
        elif self.readreverse is not None:
            self.sequencingtype = 'Paired-End'
        return self.sequencingtype

    def checkoutputpath(self):
        if self.inputistance.outputfolder is not None:
            if self.inputistance.outputfolder.endswith('/') is True:
                self.outputfolder = self.inputistance.outputfolder
            else:
                self.outputfolder = self.inputistance.outputfolder + '/'
            return self.outputfolder

    def logfilecreation(self):
        """
        Function that open new log file.
        :return: path + name file log
        """
        if os.access(self.outputfolder, os.W_OK) is True:
            self.filelog = open(self.outputfolder+self.outputid+"_mapping.log", "a")
            self.filelog.close()
            return self.outputfolder+self.outputid + "_mapping.log"
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

    # def annotesting(self, annotationfile):
    #     """
    #     check annotation
    #     :param annotationfile:
    #     :return:
    #     """
    #     self.filelog = self.openlog()
    #     if annotationfile is not None:
    #         self.annotation = annotationfile
    #         self.filelog.write(msg4 + self.annotation)
    #     else:
    #         self.filelog.write(msg12)
    #         sys.exit(1)

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
        self.filelog = self.openlog()
        if varreads is None:
            self.filelog.write(message)
            sys.exit(1)
        else:
            return varreads

    def checkfasta(self):
        """

        :return:
        """
        self.filelog = self.openlog()
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
        self.filelog = self.openlog()
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
        Log compilation
        :return:
        """
        self.filelog = self.openlog()
        self.filelog.write(datetime.datetime.now().ctime() + '\n')
        self.filelog.write(msg14)
        self.filelog.write(msg15+self.sequencingtype)
        if self.sequencingtype == 'Paired-End':
            # check dataset
            self.filelog.write(msg17 + self.checkreads(varreads=self.readforward,
                                                       message=msg2))
            self.filelog.write(msg18 + self.checkreads(varreads=self.readreverse,
                                                       message=msg3))
            self.filelog.write(msg21 + str(self.primer5forward))
            self.filelog.write(msg22 + str(self.primer3forward))
            self.filelog.write(msg23 + str(self.primer5reverse))
            self.filelog.write(msg24 + str(self.primer3reverse))
            self.filelog.write(msg47 + self.fastqcount(fastq=self.readforward,
                                                       rtype=self.readforwardtype))
            self.filelog.write(msg48 + self.fastqcount(fastq=self.readreverse,
                                                       rtype=self.readreversetype))
        elif self.sequencingtype == 'Single-End':
            self.filelog.write(msg28 + self.checkreads(varreads=self.readforward,
                                                       message=msg1))
            self.filelog.write(msg29 + self.readforwardtype)
            if self.primer5forward is None:
                pass
            else:
                self.filelog.write(msg30 + self.primer5forward)
            if self.primer3forward is None:
                pass
            else:
                self.filelog.write(msg31 + self.primer3forward)
            self.filelog.write(msg49 + self.fastqcount(fastq=self.readforward,
                                                       rtype=self.readforwardtype))
        else:
            sys.exit(1)
        self.filelog.write(msg25 + self.outputid)
        self.filelog.write(msg26 + self.fastasequence)
        # self.filelog.write(msg27 + self.annotation)
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
        self.filelog = None
        self.cloneslength = self.inputistance.minclonelength

        # check output folder
        if self.inputistance.outputfolder is not None:
            if self.inputistance.outputfolder.endswith('/') is True:
                self.outputfolder = self.inputistance.outputfolder
            else:
                self.outputfolder = self.inputistance.outputfolder + '/'
        else:
            self.filelog.write(msg9)
            sys.exit(0)
        # check output id
        if self.inputistance.outputid is not None:
            self.outputid = self.inputistance.outputid
        else:
            self.filelog.write(msg10)
            sys.exit(0)
        # check fasta sequence
        if self.inputistance.fastasequence is not None:
            self.fastasequence = self.inputistance.fastasequence
        else:
            self.filelog.write(msg11)
            sys.exit(0)
        # check annotation
        if self.inputistance.annotation is not None:
            self.annotation = self.inputistance.annotation
        else:
            self.filelog.write(msg12)
            sys.exit(0)
        # check if the name of chromosome was given by the user
        if self.inputistance.chromosomename is not None:
            self.chromosomename = self.inputistance.chromosomename
        else:
            self.filelog.write(msg87)
            sys.exit(0)
        self.out = self.outputfolder+self.outputid
        # check if log file is already created
        if self.inputistance.log is None:
            if os.access(self.outputfolder, os.W_OK) is True:
                self.filelog = open(self.outputfolder + self.outputid + ".log", "w")
                self.filelog.close()
            else:
                self.filelog.write(msg13)
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
