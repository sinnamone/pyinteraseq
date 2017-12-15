import sys
import pybedtools
import pysam
import pandas as pd
import warnings
import traceback
from output_message import *
import subprocess
from pyinteraseq_inputcheck import InputCheckDomainDefinition
import os


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
        self.filelog = self.logopen()
        try:
            # import blastnoutput
            self.df1 = pd.read_csv(blastnoutput, sep="\t", header=0,
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
            self.df3 = self.df2minus[['chr', 'end', 'start', 'seq', 'strand', 'length', 'sub']]
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
        self.filelog = self.logopen()
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
        :return: file path and name
        """
        self.filelog = self.logopen()
        self.sum = 0
        try:
            c = list(str(pysam.depth(bam)).rsplit("\n"))
            for i in range(len(c)-1):
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
        :return: flag, file path and name
        """
        self.filelog = self.logopen()
        try:
            if depb >= dept:
                self.down = float(depb) / dept
                self.bam = bamb
                self.flag = "background"
            elif depb < dept:
                self.down = float(dept) / depb
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
        :return: file path and name
        """
        self.filelog = self.logopen()
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
            self.filelog.write(msg131)
            sys.exit(1)
        else:
            self.filelog.write(msg132)
            return self.out + prefix + '.tsv'

    def intervaldata(self, bedb, bedt):
        """
        Get intervals with positive coverage
        :param bedb: bed with background coverage
        :param bedt: bed with target coverage
        :return: file path and name
        """
        self.filelog = self.logopen()
        try:
            # import background
            self.df1 = pd.read_csv(bedb, sep='\t', header=None, names=['CHR_1', 'POS_1', 'COV_1'])
            # import target
            self.df2 = pd.read_csv(bedt, sep='\t', header=None, names=['CHR_2', 'POS_2', 'COV_2'])
            # concat
            self.df3 = pd.concat([self.df1, self.df2], axis=1)
            # remove column
            self.df3.pop('CHR_2')
            # remove column
            self.df3.pop('POS_2')
            # substraction coverage
            self.df3['score'] = self.df3.COV_2 - self.df3.COV_1
            # filtering
            self.df3 = self.df3[(self.df3.COV_1 != 0) & (self.df3.score > 0)]
            # get tsv with intervals
            self.df3.to_csv(self.out + '_intervals.tsv', index=False, header=False, sep='\t')
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg133)
            sys.exit(1)
        else:
            self.filelog.write(msg134)
            return self.out + '_intervals.tsv'

    def trasposedomains(self, intervals):
        """

        :param intervals:
        :return:
        """
        self.filelog = self.logopen()
        try:
            self.df1 = pd.read_csv(intervals, sep="\t", header=None)
            self.df1[5] = self.df1[1] + 1
            self.df2 = self.df1[[0, 1, 5]]
            self.df2[1] = self.df2[1].astype(int)
            self.df2[5] = self.df2[5].astype(int)
            self.df2.to_csv(self.out + 'transposedintervals.tsv', sep="\t", header=None, index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg135)
            sys.exit(1)
        else:
            self.filelog.write(msg136)
            return self.out + 'transposedintervals.tsv'

    def filtering_domain(self, tsv):
        """
        Get intervals domains
        :param tsv: output from trasposedomains
        :return:
        """
        self.filelog = self.logopen()
        try:
            a = pybedtools.example_bedtool(tsv)
            c = a.merge(d=50)
            self.df1 = pd.read_table(c.fn, sep="\t", header=None)
            self.df1[3] = self.df1[2] - self.df1[1]
            self.df2 = self.df1.loc[self.df1[3] > 50]
            self.df3 = self.df2.loc[self.df2[3] < 1000]
            self.df3 = self.df3.rename(columns={0: 'GeneName', 1: 'Start', 2: 'End'})
            self.df3[['GeneName', 'Start', 'End']].to_csv(self.out + '_intevals_domains.txt', sep="\t", header=True, index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg137)
            sys.exit(1)
        else:
            self.filelog.write(msg138)
            return self.out + '_domains.txt'

    def cleantempfile(self):
        """
        Remove temporany files.
        :return:
        """
        templistfiles = ["background.bed", "target.bed", "background.bam", "target.bam", "target_down.bam",
                         "target.tsv", "background.tsv", "_intervals.tsv", "transposedintervals.tsv",
                         "_blastn_nohash.tab"]
        os.remove(self.outputfolder + self.namefilefasta.split('.')[0] + '.genome')
        for item in templistfiles:
            if os.path.isfile(self.out + item):
                os.remove(self.out + item)

