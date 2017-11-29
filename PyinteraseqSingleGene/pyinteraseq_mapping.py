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

    def fastq2fasta(self, fastq, nameid):
        """
        Covert Fastq in Fasta format
        :param fastq: input Fastq file
        :param nameid: output name
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
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
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
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
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
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
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
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
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
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
                                                       'fasta'))
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
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
        fnull = open(os.devnull, 'w')
        try:
            subprocess.check_call(['python', self.path_multiblastn,
                                   '--referencefasta', fasta,
                                   '--multifastasequence', multifasta,
                                   '--dbname', self.genename,
                                   '--outputfolder', self.outputfolder,
                                   '--outputid', self.outputid + suffix,
                                   '--thread', self.thread,
                                   '--outformat', outputformat,
                                   '--log', self.outputfolder + self.outputid + "_mapping.log"],
                                  stdout=fnull, stderr=fnull)
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
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
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
        self.filelog = open(self.outputfolder + self.outputid + "_mapping.log", "a")
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
            # filter in start clone
            #self.dfstart = self.dflen[(self.dflen.cstart <= 1)]
            self.dfstart = self.dflen
            # drop duplicate
            self.dfstart = self.dfstart.drop_duplicates(subset='seq', keep=False)
            if self.sequencingtype == 'Paired-End':
                # split field seq in two columns
                self.dfstart['read'], self.dfstart['seqid'] = self.dfstart['seq'].str.split(':', 2).str[0:2].str
                # split into two df read1 and read2
                self.df1 = self.dfstart[(self.dfstart['read'] == 'seq1')]
                self.df2 = self.dfstart[(self.dfstart['read'] == 'seq2')]
                self.df1['nread'] = self.df1['seq'].str.split(':', 2).str[2]
                self.df2['nread'] = self.df2['seq'].str.split(':', 2).str[2]
                # merge df
                self.dfMerge = pd.merge(self.df1, self.df2, on='nread')
                # write output
                self.dfForw = self.dfMerge[['seq_x', 'nseq_x']]
                self.dfRev = self.dfMerge[['seq_y', 'nseq_y']]
                self.dfForw = self.dfForw.rename(columns={'seq_x': 'seq', 'nseq_x': 'nseq'})
                self.dfRev = self.dfRev.rename(columns={'seq_y': 'seq', 'nseq_y': 'nseq'})
                self.dfMerge2 = self.dfForw.append(self.dfRev, ignore_index=True)
                self.dfMerge2[['seq', 'nseq']].to_csv(self.out + '_filtered_paired.tab', header=None, sep='\t',
                                                      index=False)
                self.dfMerge.to_csv(self.out + '_filtered_paired_complete.tab', header=None, sep='\t', index=False)
                return self.out + '_filtered_paired.tab'
            elif self.sequencingtype == 'Single-End':
                self.dfstart.to_csv(self.out + '_filtered_single_complete.tab', header=None, sep='\t', index=False)
                self.dfstart[['seq', 'nseq']].to_csv(self.out + '_filtered_single.tab', header=None, sep='\t',
                                                     index=False)
                return self.out + '_filtered_single.tab'
        except Warning:
            self.filelog.write('\nWarning')
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg102)
            sys.exit(1)
        else:
            self.filelog.write(msg103)

    def cleansingle(self, tempdict, sequencingtype):
        """

        :param tempdict:
        :param sequencingtype:
        :return:
        """
        if sequencingtype == "Single-End":
            tempfilelist = ["Trimmed5single", "FastaReadsForward", "TabularRenamedForward", "blastoutput",
                            "blastoutputnohash", "TabularReadsForward"]
            for i in tempfilelist:
                os.remove(tempdict[i])
            return 0
        elif sequencingtype == "Paired-End":
            tempfilelist = ["TabularRenamedForward", "FastaRenamedForward",
                            "TabularRenamedReverse", "blastoutputnohash", "FastaReadsReverse",
                            "FastaRenamedReverse", "TabularReadsReverse",
                            "TabularReadsForward", "blastoutput", "Trimmedreadconcatenated"]
            for i in tempfilelist:
                os.remove(tempdict[i])
            for i in tempdict["Trimmed5paired"]:
                os.remove(tempdict["Trimmed5paired"][i])
            return 0
        else:
            return 1

