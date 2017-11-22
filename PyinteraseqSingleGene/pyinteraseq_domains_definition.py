import sys
import pybedtools
import pysam
import os
import pandas as pd
import warnings
import traceback
from output_message import *
from Bio import SeqIO
import subprocess
from pyinteraseq_inputcheck import InputCheckDomainDefinition

# bedtools='/usr/local/bin/'
# awk='/usr/bin/awk'
# samtools='/usr/bin/samtools'


class DomainsDefinition(InputCheckDomainDefinition):

    def __init__(self, optparseinstance):
        InputCheckDomainDefinition.__init__(self, optparseinstance)
        self.df1 = None
        self.df2 = None
        self.df2plus = None
        self.df2minus = None
        self.df3 = None
        self.df4 = None
        self.prefix = None
        self.bam = None
        self.down = None
        self.flag = None
        self.strandplus = '+'
        self.strandminus = '-'
        self.depthlist = []
        self.sum = 0
        warnings.filterwarnings("ignore")

    def blastntobed(self, blastnoutput, idex):
        """
        Convert blastn output to bed format
        :param blastnoutput:
        :param idex:
        :return:
        """
        # open log
        self.filelog = open(self.outputfolder + self.outputid + "_domains_definition.log", "a")
        try:
            # import blastnoutput
            self.df1 = pd.read_csv(blastnoutput, sep="\t", header=None,
                                   names=['seq', 'chr', 'percmatch', 'length', 'mismatch', 'op', 'cstart', 'cend',
                                          'start',
                                          'end', 'evalue', 'bitscore', 'nseq'])
            # slice columns
            self.df2 = self.df1[['chr', 'start', 'end', 'seq', 'length']]
            self.df2['sub'] = self.df2['end'].sub(self.df2['start'], axis=0)
            # order dataframe
            self.df2plus = self.df2.loc[self.df2['sub'] > 0]
            self.df2minus = self.df2.loc[self.df2['sub'] < 0]
            # add strand
            self.df2plus['strand'] = self.strandplus
            self.df2minus['strand'] = self.strandminus
            # invert columns
            self.df3 = self.df2minus[['chr','end','start','seq','strand','length','sub']]
            # rename
            self.df3 = self.df3.rename(columns={'end': 'start', 'start': 'end'})
            #
            self.df4 = self.df3.append(self.df2plus, ignore_index=True)
            # sort
            self.df4.sort_values('start', inplace=True)
            # write bed file
            self.df4[['chr', 'start', 'end', 'seq', 'length', 'strand']].to_csv(self.out + idex + '.bed',
                                                                                sep="\t", header=None, index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg118)
            sys.exit(1)
        else:
            self.filelog.write(msg119)
            return self.out + idex + '.bed'

    def bedtobam(self, bedinp, genomefile, idex):
        """
        Function for convert BED to BAM file
        :param bedinp: bed file that is the output of blastntobed function
        :param genomefile: genome file that is the output of
        :param idex: string that will be appended to output name file
        :return: file path and name
        """
        self.filelog = open(self.outputfolder + self.outputid + "_domains_definition.log", "a")
        try:
            bed = pybedtools.example_bedtool(bedinp)
            bam = bed.to_bam(g=genomefile)
            bam.moveto(self.out + idex + '.bam')
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg125)
            sys.exit(1)
        else:
            self.filelog.write(msg126)
            return self.out + idex + '.bam'

    def averagedepth(self, bam, lung):
        """
        Computing average depth BAM file
        :param bam:
        :param lung:
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_domains_definition.log", "a")
        try:
            c = list(str(pysam.depth(bam)).rsplit("\n"))
            for i in range(len(c)):
                self.depthlist.append(c[i][0:])
            for g in range(len(self.depthlist) - 1):
                self.sum += int(self.depthlist[g][-1])
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg127)
            sys.exit(1)
        else:
            self.filelog.write(msg128)
            return float(self.sum) / lung

    def downsampling(self, depb, dept, bamb, bamt):
        """

        :param depb:
        :param dept:
        :param bamb:
        :param bamt:
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_domains_definition.log", "a")
        try:
            if depb >= dept:
                self.down = float(dept) / depb
                self.bam = bamb
                self.flag = "background"
            elif dept >= depb:
                self.down = float(depb) / dept
                self.bam = bamt
                self.flag = "target"
            self.prefix = self.bam.rsplit("/")[-1][:-len(".bam")]
            with open(self.outputfolder + self.prefix + '_down.bam', 'w') as a:
                subprocess.check_call(['samtools', 'view', '-bs', str(self.down), self.bam], stdout=a)
        except subprocess.CalledProcessError:
            self.filelog.write(msg129)
            sys.exit(1)
        else:
            self.filelog.write(msg130)
            return self.flag, self.outputfolder + self.prefix + '_down.bam'

    def depthcount(self, bam, prefix):
        """
        Coverage computation
        :param bam: bam file
        :param prefix: id to append to output
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + "_domains_definition.log", "a")
        try:
            # read bam file  output all positions
            self.depthlist = pysam.depth("-a", bam).split("\n")
            # import list in dataframe
            self.df1 = pd.DataFrame(data=[self.depthlist]).transpose()
            # split columns into 3 colums
            self.df2 = self.df1[0].apply(lambda x: pd.Series(str(x).split('\t')))
            # save
            self.df2.to_csv(self.out + prefix + '.tsv', sep="\t", header=None, index=False)
        except traceback:
            self.filelog.write(msg129)
            sys.exit(1)
        else:
            self.filelog.write(msg130)
            return self.out + prefix + '.tsv'
