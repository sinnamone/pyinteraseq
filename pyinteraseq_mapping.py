import subprocess
from output_message import *
from Bio import SeqIO
import os
import pandas as pd
import sys
from pyinteraseq_inputcheck import InputCheck
from multiprocessing import Pool


class BlastNlucleotide(InputCheck):

    def __init__(self, optparseinstance):
        InputCheck.__init__(self, optparseinstance)
        self.out_lines = []
        self.temp_line = ''
        self.df1 = None
        self.seqix = 1
        self.id = '1'
        self.id = [(int(x) - 1) for x in self.id]
        self.header = None
        self.dbid = os.path.basename(self.fastasequence.split('/')[-1])
        self.path_multiblastn = os.path.dirname(os.path.realpath(__file__)) + '/pyinteraseq_multblastn.py'
        self.pool = Pool(processes=int(self.thread))

    def fastq2fasta(self, fastq, nameid):
        """
        Covert Fastq in Fasta format
        :param fastq: input Fastq file
        :param nameid: output name
        :return:
        """
        if nameid == "forward":
            SeqIO.convert(fastq, 'fastq', self.out + '1.fasta', 'fasta')
            return self.out + '1.fasta'
        elif nameid == "reverse":
            SeqIO.convert(fastq, 'fastq', self.out + '1.fasta', 'fasta')
            return self.out + '1.fasta'

    def seqrename(self, tabular, readirection):
        if readirection == "forward":
            self.df1 = pd.read_csv(tabular, header=None, sep='\t')
            self.df1['seq_id'] = self.df1.apply(lambda x: "seq1:1:" + str(x.name), axis=1)
            self.df1[['seq_id', 1]].to_csv(self.out + '_1_newid.tab', header=None, sep='\t', index=False)
            return self.out + '_1_newid.tab'
        else:
            self.df1 = pd.read_csv(tabular, header=None, sep='\t')
            self.df1['seq_id'] = self.df1.apply(lambda x: "seq2:2:" + str(x.name), axis=1)
            self.df1[['seq_id', 1]].to_csv(self.out + '_2_newid.tab', header=None, sep='\t', index=False)
            return self.out + '_2_newid.tab'

    def tab2fasta(self, tabular, prefixoutput):
        with open(tabular, 'r') as f:
            with open(self.out + prefixoutput + '.fasta', 'w') as f_out:
                for line in f:
                    line = line.strip().split('\t')
                    self.header = '>' + '_'.join([line[i] for i in self.id])
                    f_out.write(self.header + '\n')
                    f_out.write(line[self.seqix] + '\n')
            f_out.close()
        return self.out + prefixoutput + '.fasta'

    def concatenateforrev(self, readlist):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        with open(self.out + '_con.fasta', 'w') as outfile:
            for fname in readlist:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
        outfile.close()
        # self.filelog.write(msg56 + self.fastqcount(self.out + '_con.fasta',
        #                                            self.readforwardtype))
        return self.out + '_con.fasta'

    def fasta2tabular(self, imp, prefix):
        SeqIO.convert(imp, 'fasta', self.out + prefix + '.tab', 'tab')
        return self.out + prefix + '.tab'

    def callmultiblastn(self, fasta, multifasta, outputformat, suffix):
        print fasta
        print multifasta
        print outputformat
        print suffix
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            subprocess.check_call(['python', self.path_multiblastn,
                                   '--referencefasta', fasta,
                                   '--multifastasequence', multifasta,
                                   '--dbname', self.chromosomename,
                                   '--outputfolder', self.outputfolder,
                                   '--outputid', self.outputid + suffix,
                                   '--thread', self.thread,
                                   '--outformat', outputformat])
        except subprocess.CalledProcessError:
            self.filelog.write(msg60)
            sys.exit(0)
        else:
            self.filelog.write(msg61)
            return self.out + suffix
