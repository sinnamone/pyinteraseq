from output_message import *
import sys
import os
import datetime
import subprocess
from Bio import SeqIO
import traceback


class InputCheck(object):

    def __init__(self, optparseinstance):
        # import instance with all input flag
        self.inputistance = optparseinstance
        # put all under this
        self.thread = self.inputistance.thread
        self.count = 0
        self.filelog = None
        self.cutadapt = ''
        self.cutadaptversion = ""
        self.cloneslength = self.inputistance.minclonelength
        self.readforward = self.inputistance.readforward
        self.readreverse = self.inputistance.readreverse
        self.sequencingtype = None
        # self.readforwardtype = self.inputistance.readforwardtype
        # self.readreversetype = self.inputistance.readreversetype
        self.readforwardtype = None
        self.readreversetype = None
        self.primer5forward = self.inputistance.primer5forward
        self.primer3forward = self.inputistance.primer3forward
        self.primer5reverse = self.inputistance.primer5reverse
        self.primer3reverse = self.inputistance.primer3reverse
        self.fastasequence = self.inputistance.fastasequence
        self.genename = self.inputistance.genename
        # self.outputfolder+self.outputid
        #self.awk = None
        self.pick_otus = None
        self.pick_rep_set = None
        self.samtools = None
        self.outputfolder = None
        self.outputid = self.inputistance.outputid
        self.out = None
        #
        if self.inputistance.outputfolder is not None:
            if self.inputistance.outputfolder.endswith('/') is True:
                self.outputfolder = self.inputistance.outputfolder
                self.out = self.outputfolder + self.outputid
            else:
                self.outputfolder = self.inputistance.outputfolder + '/'
                self.out = self.outputfolder + self.outputid
        else:
            sys.exit(1)
        #
        if self.inputistance.outputid is not None:
            self.outputid = self.inputistance.outputid

    def logfilecreation(self):
        """
        Function that open new log file.
        :return: path + name file log
        """
        if os.access(self.outputfolder, os.W_OK) is True:
            self.filelog = open(self.outputfolder+self.outputid+"_trimming.log", "w")
            self.filelog.close()
            return self.outputfolder+self.outputid+"_trimming.log"
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
        self.filelog = open(self.outputfolder + self.outputid + "_trimming.log", "a")
        if varreads is None:
            self.filelog.write(message)
            sys.exit(1)
        else:
            return varreads

    def checkreadtype(self, varformattpye):
        """
        Check if read input is fasta or fastq
        :param varformattpye:
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_trimming.log", "a")
        # Check type of dataset
        if varformattpye is None:
            self.filelog.write(msg4)
            sys.exit(1)
        else:
            return varformattpye

    def checksampletype(self, varsampletype, outputmessage):
        """

        :param varsampletype:
        :param outputmessage:
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_trimming.log", "a")
        # Check for Dataset typr
        if varsampletype is None:
            self.filelog.write(outputmessage)
            sys.exit(1)
        else:
            return varsampletype

    def checkprimers(self, varprimers, message):
        """

        :param varprimers:
        :param message:
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_trimming.log", "a")
        # Check primers
        if varprimers is None:
            sys.stdout.write(message)
            sys.exit(1)
        else:
            return varprimers

    def checkfasta(self):
        """

        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_trimming.log", "a")
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
        self.filelog = open(self.outputfolder + self.outputid + "_trimming.log", "a")
        # check if the name of chromosome was given by the user
        if self.inputistance.genename is None:
            sys.stdout.write(msg87)
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

    def samtoolscheck(self):
        """
        Function that verify the installation of samtools
        :return: samtools path
        """
        try:
            self.samtools = subprocess.check_output(['which', 'samtools']).split('\n')[0]
        except subprocess.CalledProcessError:
            self.filelog.write(msg69)
            sys.exit(1)
        else:
            return self.samtools

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

        self.filelog = open(self.outputfolder + self.outputid + "_trimming.log", "a")
        self.filelog.write(datetime.datetime.now().ctime() + '\n')
        self.filelog.write(msg14)
        self.filelog.write(msg15 + self.checkreversereads(readreverse=self.readreverse))
        if self.sequencingtype == 'Paired-End':
            # check dataset
            self.filelog.write(msg17 + self.checkreads(varreads=self.readforward,
                                                       message=msg2))
            self.filelog.write(msg18 + self.checkreads(varreads=self.readreverse,
                                                       message=msg3))
            # check type
            self.readforwardtype = self.fastatesting(self.readforward)
            self.readreversetype = self.fastatesting(self.readreverse)
            self.filelog.write(msg19 + self.checkreadtype(varformattpye=self.readforwardtype))
            self.filelog.write(msg20 + self.checkreadtype(varformattpye=self.readreversetype))
            # check primers
            self.filelog.write(msg21 + self.checkprimers(varprimers=self.primer5forward,
                                                         message=msg6))
            self.filelog.write(msg22 + self.checkprimers(varprimers=self.primer3forward,
                                                         message=msg5))
            self.filelog.write(msg23 + self.checkprimers(varprimers=self.primer5reverse,
                                                         message=msg8))
            self.filelog.write(msg24 + self.checkprimers(varprimers=self.primer3reverse,
                                                         message=msg7))
            self.filelog.write(msg47 + self.fastqcount(fastq=self.readforward,
                                                       rtype=self.readforwardtype))
            self.filelog.write(msg48 + self.fastqcount(fastq=self.readreverse,
                                                       rtype=self.readreversetype))
        else:
            self.filelog.write(msg28 + self.checkreads(varreads=self.readforward,
                                                       message=msg1))
            self.filelog.write(msg29 + self.checkreadtype(varformattpye=self.readforwardtype))
            self.filelog.write(msg30 + self.checkprimers(varprimers=self.primer5forward,
                                                         message=msg6))
            self.filelog.write(msg31 + self.checkprimers(varprimers=self.primer3forward,
                                                         message=msg5))
            self.filelog.write(msg49 + self.fastqcount(self.checkprimers(varprimers=self.primer5forward,
                                                                         message=msg6),
                                                       self.checkprimers(varprimers=self.primer3forward,
                                                                         message=msg5)
                                                       ))
        self.filelog.write(msg25 + self.outputid)
        self.filelog.write(msg26 + self.checkfasta())
        self.filelog.write(msg27 + self.checkgenename())
        self.filelog.write(msg33 + subprocess.check_output(['cutadapt', '--version']))
        self.filelog.write(subprocess.check_output(['pick_otus.py', '--version']))
        self.filelog.write(subprocess.check_output(['pick_rep_set.py', '--version']))
        self.filelog.write(msg0)
        self.filelog.close()
        return True


class InputCheckMapping(object):

    def __init__(self, optparseinstance):
        self.inputistance = optparseinstance
        self.readforwardtrimmed = self.inputistance.readforwardtrimmed
        self.readreversetrimmed = self.inputistance.readreversetrimmed
        self.readforwardtrimmedtype = self.inputistance.readforwardtrimmedtype
        self.readreversetrimmedtype = self.inputistance.readreversetrimmedtype
        self.thread = self.inputistance.thread
        self.count = 0
        self.filelog = None
        self.cloneslength = self.inputistance.minclonelength
        self.sequencingtype = None
        self.fastasequence = self.inputistance.fastasequence
        self.genename = self.inputistance.genename
        # self.outputfolder+self.outputid
        self.outputfolder = None
        self.outputid = self.inputistance.outputid
        self.out = None
        self.opengap = self.inputistance.opengap
        self.mismatch = self.inputistance.mismatch
        #
        if self.inputistance.outputfolder is not None:
            if self.inputistance.outputfolder.endswith('/') is True:
                self.outputfolder = self.inputistance.outputfolder
                self.out = self.outputfolder + self.outputid
            else:
                self.outputfolder = self.inputistance.outputfolder + '/'
                self.out = self.outputfolder + self.outputid
        else:
            sys.exit(1)
        #
        if self.inputistance.outputid is not None:
            self.outputid = self.inputistance.outputid

    def inputinformationappen(self):
        """
        Log compilation. Call all the previous functions
        :return:
        """

        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "w")
        self.filelog.write(datetime.datetime.now().ctime() + '\n')
        self.filelog.write(msg14)
        self.filelog.write(msg15 + self.checkreversereads(readreverse=self.readreversetrimmed))
        if self.sequencingtype == 'Paired-End':
            # check dataset
            self.filelog.write(msg17 + self.checkreads(varreads=self.readforwardtrimmed,
                                                       message=msg2))
            self.filelog.write(msg18 + self.checkreads(varreads=self.readreversetrimmed,
                                                       message=msg3))
            # check type
            self.filelog.write(msg19 + self.checkreadtype(varformattpye=self.readforwardtrimmedtype))
            self.filelog.write(msg20 + self.checkreadtype(varformattpye=self.readreversetrimmedtype))
            self.filelog.write(msg47 + self.fastqcount(fastq=self.readforwardtrimmed,
                                                       rtype=self.readforwardtrimmedtype))
            self.filelog.write(msg48 + self.fastqcount(fastq=self.readreversetrimmed,
                                                       rtype=self.readreversetrimmedtype))
        else:
            self.filelog.write(msg28 + self.checkreads(varreads=self.readforwardtrimmed,
                                                       message=msg1))
            self.filelog.write(msg29 + self.checkreadtype(varformattpye=self.readforwardtrimmedtype))
            self.filelog.write(msg49 + self.fastqcount(fastq=self.readforwardtrimmed,
                                                       rtype=self.readforwardtrimmedtype))
        self.filelog.write(msg25 + self.outputid)
        self.filelog.write(msg26 + self.checkfasta())
        self.filelog.write(msg27 + self.checkgenename())
        self.filelog.write(msg33 + subprocess.check_output(['cutadapt', '--version']))
        self.filelog.write(subprocess.check_output(['pick_otus.py', '--version']))
        self.filelog.write(subprocess.check_output(['pick_rep_set.py', '--version']))
        self.filelog.write(msg0)
        self.filelog.close()
        return True

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
            self.readreversetrimmed = self.inputistance.readreversetrimmed
            self.sequencingtype = 'Paired-End'
            return self.sequencingtype

    def checkreads(self, varreads, message):
        """
        Chech if forward or reverse reads exist
        :param varreads:
        :param message:
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
        if varreads is None:
            self.filelog.write(message)
            sys.exit(1)
        else:
            return varreads

    def checkreadtype(self, varformattpye):
        """
        Check if read input is fasta or fastq
        :param varformattpye:
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
        # Check type of dataset
        if varformattpye is None:
            self.filelog.write(msg4)
            sys.exit(1)
        else:
            return varformattpye

    def checksampletype(self, varsampletype, outputmessage):
        """

        :param varsampletype:
        :param outputmessage:
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
        # Check for Dataset typr
        if varsampletype is None:
            self.filelog.write(outputmessage)
            sys.exit(1)
        else:
            return varsampletype

    def checkfasta(self):
        """

        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
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
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
        # check if the name of chromosome was given by the user
        if self.inputistance.genename is None:
            sys.stdout.write(msg87)
            sys.exit(1)
        else:
            return self.genename

    def logcheck(self):
        """

        :return:
        """
        if self.inputistance.log is None:
            if os.access(self.outputfolder, os.W_OK) is True:
                self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "w")
                self.filelog.close()
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


class InputCheckDomainDefinition(object):

    def __init__(self, optparseinstance):
        # import instance with all input flag
        self.inputistance = optparseinstance
        self.count = 0
        self.filelog = None
        self.genome = None
        self.backgroundmappingoutput = self.inputistance.backgroundmappingoutput
        self.targetmappingoutput = self.inputistance.targetmappingoutput
        # check output folder
        if self.inputistance.outputfolder is not None:
            if self.inputistance.outputfolder.endswith('/') is True:
                self.outputfolder = self.inputistance.outputfolder
            else:
                self.outputfolder = self.inputistance.outputfolder + '/'
        else:
            sys.stdout.write(msg9)
            sys.exit(1)
        # check output id
        if self.inputistance.outputid is not None:
            self.outputid = self.inputistance.outputid
            self.out = self.outputfolder + self.outputid
        else:
            sys.stdout.write(msg10)
            sys.exit(1)
        # check fasta sequence
        if self.inputistance.fastasequence is not None:
            self.fastasequence = self.inputistance.fastasequence
        else:
            sys.stdout.write(msg11)
            sys.exit(1)
        if self.inputistance.genename is None:
            sys.stdout.write(msg11)
            sys.exit(1)
        else:
            self.genename = self.inputistance.genename

    def inputinformationappen(self):
        """
        Log compilation. Call all the previous functions
        :return:
        """

        self.filelog = open(self.outputfolder + self.outputid + "_domains_definition.log", "w")
        self.filelog.write(datetime.datetime.now().ctime() + '\n')
        self.filelog.write(msg14)
        self.filelog.write(msg114 + self.backgroundmappingoutput)
        self.filelog.write(msg115 + self.targetmappingoutput)
        self.filelog.write(msg116 + self.outputid)
        self.filelog.write(msg117 + self.outputid)
        self.filelog.close()
        return True

    def contafasta(self, ref):
        """
        Count gene lenght
        :param ref:
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_domains_definition.log", "a")
        try:
            for seq_record in SeqIO.parse(ref, "fasta"):
                leng = len(seq_record)
                return leng
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg121)
            sys.exit(1)
        else:
            self.filelog.write(msg122)

    def genomefile(self):
        """
        Function for the creation of genome file
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_domains_definition.log", "a")
        try:
            ln = self.contafasta(ref=self.fastasequence)
            self.genome = open(self.outputfolder + self.genename + ".genome", "w")
            self.genome.write(self.genename + '\t' + '1\t' + str(ln))
            self.genome.close()
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg123)
            sys.exit(1)
        else:
            self.filelog.write(msg124)
            return self.outputfolder + self.genename + ".genome"
