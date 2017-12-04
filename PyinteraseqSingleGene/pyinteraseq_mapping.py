import subprocess
from output_message import *
from Bio import SeqIO
import os
import pandas as pd
import sys
from pyinteraseq_inputcheck import InputCheck
from multiprocessing import Pool
import traceback
import warnings


class BlastNlucleotide(InputCheck):

    def __init__(self, optparseinstance):
        InputCheck.__init__(self, optparseinstance)
        warnings.filterwarnings("ignore")
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
        self.df = None
        self.df2 = None
        self.dfOp = None
        self.dfMM = None
        self.dflen = None
        self.dfstart = None
        self.dfMerge = None
        self.dfForw = None
        self.dfRev = None
        self.dfMerge2 = None
        self.dbname = self.outputfolder + os.path.basename(self.fastasequence.split('/')[-1]).split('.')[0]



    def fastq2fasta(self, fastq, nameid):
        """
        Covert Fastq in Fasta format
        :param fastq: input Fastq file
        :param nameid: output name
        :return:
        """
        self.filelog = self.logopen()
        try:
            SeqIO.convert(fastq, 'fastq', self.out + nameid + '.fasta', 'fasta')
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg94)
            sys.exit(1)
        else:
            self.filelog.write(msg95)
            return self.out + nameid + '.fasta'

    def fasta2tabular(self, imp, prefix):
        """
        Function for conversion Fasta file in tabular format
        :param imp: input file in fasta format
        :param prefix: prefix to add at converted file
        :return:
        """
        self.filelog = self.logopen()
        try:
            SeqIO.convert(imp, 'fasta', self.out + prefix + '.tab', 'tab')
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg96)
            sys.exit(1)
        else:
            self.filelog.write(msg97)
            return self.out + prefix + '.tab'

    def seqrename(self, tabular, readirection):
        """
        Rename Id fasta, first sequence will be 0 last is seq count
        :param tabular: File fasta in tabular format
        :param readirection: "forward" or "reverse"
        :return: path + name + _1_newid.tab
        """
        self.filelog = self.logopen()
        try:
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
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg98)
            sys.exit(1)
        finally:
            self.filelog.write(msg99)

    def tab2fasta(self, tabular, prefixoutput):
        """
        Convert tabular fasta file (2 columns) first ID second nucleotide sequence in fasta file format
        :param tabular: File fasta in tabular format
        :param prefixoutput: prefix to append
        :return: path + idanalysis + prefix + .fasta
        """
        self.filelog = self.logopen()
        try:
            with open(tabular, 'r') as f:
                with open(self.out + prefixoutput + '.fasta', 'w') as f_out:
                    for line in f:
                        line = line.strip().split('\t')
                        self.header = '>' + '_'.join([line[i] for i in self.id])
                        f_out.write(self.header + '\n')
                        f_out.write(line[self.seqix] + '\n')
                f_out.close()
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg100)
            sys.exit(1)
        else:
            self.filelog.write(msg101)
            return self.out + prefixoutput + '.fasta'

    def concatenateforrev(self, readlist):
        """
        Merge forward and reverse fasta file
        :param readlist: list with files to append
        :return: path + idanalysis + _con.fasta
        """
        self.filelog = self.logopen()
        try:
            with open(self.out + '_con.fasta', 'w') as outfile:
                for fname in readlist:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
            outfile.close()
        except StandardError:
            self.filelog.write(msg88)
            sys.exit(1)
        else:
            self.filelog.write(msg56 + self.fastqcount(self.out + '_con.fasta',
                                                       'fasta') + '\n')
            return self.out + '_con.fasta'

    def callmultiblastn(self, fasta, multifasta, outputformat, suffix):
        """
        Function to call multithread blastn,
        :param fasta: Reference or fasta that will be use for makeblastdb
        :param multifasta: Input multifasta file
        :param outputformat: Blast output format
        :param suffix: String added to outputfile
        :return: blastn output
        """
        self.filelog = self.logopen()
        fnull = open(os.devnull, 'w')
        try:
            subprocess.check_call(['python', self.path_multiblastn,
                                   '--referencefasta', fasta,
                                   '--multifastasequence', multifasta,
                                   '--dbname', self.genename,
                                   '--outputfolder', self.outputfolder,
                                   '--outputid', self.outputid + suffix,
                                   '--thread', self.thread,
                                   '--outformat', outputformat])
        except subprocess.CalledProcessError:
            self.filelog.write(msg60)
            sys.exit(1)
        else:
            self.filelog.write(msg61)
            self.filelog.close()
            return self.out + suffix

    def hashclean(self, blastnout, prefix):
        """
        Function to clean hash in blastn 7 format output
        :param blastnout: Blastn output
        :param prefix: prefix add to output file
        :return: path + prefix + '.tab' of new file
        """
        self.filelog = self.logopen()
        try:
            with open(blastnout) as oldfile, open(self.out + prefix + '.tab', 'w') as newfile:
                for line in oldfile:
                    if not line.startswith('#'):
                        newfile.write(line)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg92)
            sys.exit(0)
        else:
            self.filelog.write(msg93)
            return self.out + prefix + '.tab'

    def blastnfiltering(self, blastnout):
        """
        Function that takes as input blastn no hash created with function hashclean and filter out , 1 opengap,
         5 %of mismatch
        :param blastnout:
        :return:
        """
        self.filelog = self.logopen()
        try:
            self.df = pd.read_csv(blastnout, sep='\t', header=None,
                                  names=['seq', 'chr', 'percmatch', 'length', 'mismatch', 'op', 'cstart', 'cend',
                                         'start', 'end', 'evalue', 'bitscore', 'nseq'])
            self.filelog.write(msg62 + str(len(self.df)))
            # filter open gap
            self.dfOp = self.df[(self.df.op <= int(self.opengap))]
            self.filelog.write(msg63 + str(len(self.dfOp)))
            # trasform mismatch in percentage of mismatch using clone length
            self.dfOp['pmismatch'] = (self.dfOp.mismatch.div(self.dfOp.length).mul(100))
            # trasform into numeric field
            self.dfOp[['pmismatch']].apply(pd.to_numeric)
            # filter on percentage of mismatch
            self.dfMM = self.dfOp[(self.dfOp['pmismatch'] < float(self.mismatch))]
            self.filelog.write(msg64 + str(len(self.dfMM)))
            # lenght filtering
            self.dflen = self.dfMM[(self.dfMM.length >= int(self.cloneslength))]
            self.filelog.write(msg65 + str(len(self.dflen)))
            # drop duplicate
            self.dfstart = self.dflen.drop_duplicates(subset='seq', keep=False)
            if self.sequencingtype == 'Paired-End':
                # split field seq in two columns
                self.dfstart['read'], self.dfstart['seqid'] = self.dfstart['seq'].str.split(':', 2).str[0:2].str
                # split into two df read1 and read2
                self.df1 = self.dfstart[(self.dfstart['read'] == 'seq1')]
                self.df2 = self.dfstart[(self.dfstart['read'] == 'seq2')]
                self.df1['nread'] = self.df1['seq'].str.split(':', 2).str[2]
                self.df2['nread'] = self.df2['seq'].str.split(':', 2).str[2]
                # merge df
                #self.dfMerge = pd.merge(self.df1, self.df2, on='nread')
                # write output
                # self.dfForw = self.dfMerge[['seq_x', 'nseq_x']]
                # self.dfRev = self.dfMerge[['seq_y', 'nseq_y']]
                # self.dfForw = self.dfForw.rename(columns={'seq_x': 'seq', 'nseq_x': 'nseq'})
                # self.dfRev = self.dfRev.rename(columns={'seq_y': 'seq', 'nseq_y': 'nseq'})
                self.dfMerge2 = self.df1.append(self.df2, ignore_index=True)
                self.dfMerge2[[u'seq', u'chr', u'percmatch', u'length', u'pmismatch', u'op', u'cstart',
                               u'cend', u'start', u'end', u'evalue', u'nread', u'nseq']].to_csv(self.out + '_output_mapping_step.tab', header=None, sep='\t',
                                                      index=False)
                #self.dfMerge.to_csv(self.out + '_filtered_paired_complete.tab', header=None, sep='\t', index=False)
                return self.out + '_filtered_paired.tab'
            elif self.sequencingtype == 'Single-End':
                # in case will need fasta
                #self.dfstart.to_csv(self.out + '_filtered_single_complete.tab', header=None, sep='\t', index=False)
                self.dfstart[[u'seq', u'chr', u'percmatch', u'length', u'pmismatch', u'op', u'cstart',
                              u'cend', u'start', u'end', u'evalue', u'bitscore', u'nseq']].to_csv(
                    self.out + '_output_mapping_step.tab', header=None, sep='\t', index=False)
                return self.out + '_output_mapping_step.tab'
        except Warning:
            self.filelog.write('\nWarning')
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg102)
            sys.exit(1)
        else:
            self.filelog.write(msg103)

    def cleantempfile(self):
        """
        Remove temporany files.
        :return:
        """
        db1 = str(self.dbname + ".nsq")
        db2 = str(self.dbname + ".nin")
        db3 = str(self.dbname + ".nhr")
        os.remove(db1)
        os.remove(db2)
        os.remove(db3)
        # "_filtered_paired_complete.tab"
        templistfilesingle = ["_forward.fastq", "_read1.fastq", "_forward.tab", "_forward.fasta", "_filtered_single.tab",
                              "_blastn_nohash.tab", "_blastn.txt", "_1_newid.tab"]
        templistfilepaired = ["_forward.fastq", "_reverse.fastq", "_read1.fastq", "_read2.fastq", "_forward.tab",
                              "_reverse.tab", "_1_newid.tab", "_2_newid.tab", "_forward.fasta", "_reverse.fasta",
                              "_con.fasta", "_blastn.txt", "_filtered_paired.tab", "_blastn_nohash.tab"]
        if self.sequencingtype in "Single-End":
            if self.readforwardtype in "fastq":
                for item in templistfilesingle:
                    if os.path.isfile(self.out + item):
                        os.remove(self.out + item)
            elif self.readforwardtype in "fasta":
                templistfilesingle = templistfilesingle[2:]
                templistfilesingle.append("_read1.fasta")
                for item in templistfilesingle:
                    if os.path.isfile(self.out + item):
                        os.remove(self.out + item)
        elif self.sequencingtype in "Paired-End":
            if self.readforwardtype in "fastq":
                for item in templistfilepaired:
                    if os.path.isfile(self.out + item):
                        os.remove(self.out + item)
            elif self.readforwardtype in "fasta":
                templistfilepaired = templistfilepaired[2:]
                templistfilepaired.append("_read1.fasta")
                templistfilepaired.append("_read2.fasta")
                for item in templistfilepaired:
                    if os.path.isfile(self.out + item):
                        os.remove(self.out + item)
