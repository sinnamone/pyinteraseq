import os
from os import listdir
import subprocess
import sys
import time
import pandas as pd
import numpy as np
from output_message import *
from pyinteraseq_inputcheck import InputCheckDomainDefinition
import pybedtools
import traceback
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

class DomainsDefinition(InputCheckDomainDefinition):

    def __init__(self, optparseinstance):
        InputCheckDomainDefinition.__init__(self, optparseinstance)
        self.filelog = self.logopen()

    def depthcoverage(self):
        """
        bedtoools genomecov for each base covered
        :return:
        """

        try:
            with open(self.out + '_DepthCoverageBed.txt', 'w') as b:
                subprocess.check_call(['bedtools', 'genomecov', '-dz', '-ibam', self.bamfile, '-g', self.genome],
                                      stdout=b)
        except subprocess.CalledProcessError:
            self.filelog.write(msg89)
            sys.exit(1)
        except traceback:
            self.filelog.write(msg89)
            sys.exit(1)
        else:
            self.filelog.write(msg90)
            return self.out + '_DepthCoverageBed.txt'

    def breadthcoverage(self):
        """
        Bedtools genomecov for read count
        :return:
        """
        try:
            with open(self.out + '_BreadthCoverageBed.txt', 'w') as b:
                subprocess.check_call(['bedtools', 'genomecov','-ibam', self.bamfile, '-g', self.genome],stdout=b)
        except subprocess.CalledProcessError:
            self.filelog.write(msg91)
            sys.exit(1)
        except traceback:
            self.filelog.write(msg91)
            sys.exit(1)
        else:
            self.filelog.write(msg92)
            return self.out + '_BreadthCoverageBed.txt'

    def bam2tabular(self):
        """
        Convert BAM file to Bed
        :return:
        """
        try:
            bedfile = pybedtools.BedTool(self.bamfile).bam_to_bed()
            # bedfile = pybedtools.BedTool(self.bamfile).bam_to_bed().saveas(self.out + '.bed')
            bedfiledf = pd.read_table(bedfile.fn, header=None,sep="\t")
            annot = pd.read_csv(self.annotation, header=None, sep="\t")
            df3 = pd.merge(bedfiledf, annot, left_on=0,right_on=3)
            df3["start"] = df3.apply(lambda x: x["1_y"] + x["1_x"], axis=1)
            df3["end"] = df3.apply(lambda x: x["1_y"] + x["2_x"], axis=1)
            df3[["0_y", "start", "end", "0_x", "4_x", "5_x"]].to_csv(self.out + '.bed', sep="\t", header=None, index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg94)
            sys.exit(1)
        else:
            self.filelog.write(msg93)
            return self.out + '.bed'

    def groupbyreads(self, convertedbed):
        try:
            df = pd.read_csv(convertedbed, header=None, sep='\t',names=['chr', 'start', 'stop', 'id', 'score', 'strand'])
            df1 = df.groupby(["chr"], as_index=False)["id"].count()
            df1.to_csv(self.out + '_readgroupedbychr.txt', header=None, sep='\t', index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg95)
            sys.exit(1)
        else:
            self.filelog.write(msg96)
            return self.out + '_readgroupedbychr.txt'

    def coordinatesfix(self, convertedbam):
        """
        Fix start and end and prepare tab to be annotated
        :param bedparsed:
        :return:
        """
        try:
            self.df1 = pd.read_csv(convertedbam, sep="\t", header=None,
                                   names=["ID", "start_clone", "end_clone", "id_clone", "score", "strand"])
            # print self.df1.head()

            # self.aux = self.df1["ID"].apply(lambda x: x.split('_'))
            # self.df1["start"] = self.aux.apply(lambda x: x[2]).astype(int)
            # self.df1["end"] = self.aux.apply(lambda x: x[3]).astype(int)
            # self.df1["new_start"] = self.df1.apply(lambda row: row.start + row.start_clone, axis=1)
            # self.df1["new_end"] = self.df1.apply(lambda row: row.new_start + row.end_clone, axis=1)
            # self.df1[["ID", "new_start", "new_end", "id_clone", "score", "strand"]].to_csv(self.out + "_bedfixed.tab", header=None, sep="\t", index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg74)
            sys.exit(1)
        else:
            self.filelog.write(msg73)
            return self.out + "_bedfixed.tab"

    def bedtoolscoverage(self, outputcoordinatesfix):
        try:
            bed = pybedtools.BedTool(outputcoordinatesfix)
            annotation = pybedtools.BedTool(self.annotation)
            cov = annotation.coverage(bed, counts=True).saveas(self.out + '_readcov.bed')
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg98)
            sys.exit(1)
        else:
            self.filelog.write(msg97)
            return self.out + '_readcov.bed'

    def groupbydepth(self, outputbreadthcoverage):
        try:
            df = pd.read_csv(outputbreadthcoverage, header=None, sep='\t', names=["chr", "depth", "numbasdepth", "lenght", "perccov"])
            df1 = df.groupby(["chr"],sort=False)["depth"].max()
            df1.to_csv(self.out + "_maxdepth.txt", header=None, sep='\t', index=True)
        except:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg99)
            sys.exit(1)
        else:
            self.filelog.write(msg100)
            return self.out + "_maxdepth.txt"

    def percentile(self, file):
        """
        Filtering domains with low coverage
        :param file:
        :return:
        """
        try:
            # dataframe creation
            dfA = pd.read_csv(file, index_col=False, header=None, sep='\t')
            # labeling columns
            dfA.columns = ['a1', 'a2', 'a3']
            # list with values of depth
            p = dfA.a3
            # percentile function, value int perc is the threshold
            v = np.percentile(p, self.threshold)
            dfD = dfA.loc[lambda df: df.a3 > v, :]
            dfD.to_csv(self.out + '_percentilefiltered.bed', header=None,
                       sep='\t', index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg101)
            sys.exit(1)
        else:
            self.filelog.write(msg102)
            return self.out + '_percentilefiltered.bed'

    def startenddefinition(self, file):
        """
        Setepu start and end domain
        :param file:
        :return:
        """
        try:
            df = pd.read_csv(file, index_col=False, header=None, sep='\t')
            df.columns = ['a', 'b', 'c']
            df1 = df.groupby(['a'], as_index=False)['b'].min()
            df2 = df.groupby(["a"], as_index=False)["b"].max()
            df3 = df.groupby(["a"], as_index=False)["c"].mean()
            dfA = pd.merge(df1, df2, on='a')
            dfB = pd.merge(dfA, df3, on='a')
            dfB.to_csv(self.out + '_raw_domains.txt', header=None, sep='\t', index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg103)
            sys.exit(1)
        else:
            self.filelog.write(msg104)
            return self.out + '_raw_domains.txt'

    def domainparsing(self, outputstartenddefinition):
        """

        :param outputstartenddefinition:
        :return:
        """
        try:
            df = pd.read_csv(outputstartenddefinition, sep="\t", header=None,
                             names=["ID", "start_clone", "end_clone", "ave_depth"])
            annot = pd.read_csv(self.annotation, header=None, sep="\t")
            aux = pd.merge(df, annot, left_on="ID",right_on=3)
            aux["start"] = aux.apply(lambda x: x[1] + x["start_clone"], axis=1)
            aux["end"] = aux.apply(lambda x: x[1] + x["end_clone"], axis=1)
            aux["lenght"] = aux["end"] - aux["start"]
            # aux = df["ID"].apply(lambda x: x.split('_'))
            # df["geneid"] = aux.apply(lambda x: x[0])
            # df["chr"] = aux.apply(lambda x: x[1])
            # df["start"] = aux.apply(lambda x: x[2]).astype(int)
            # df["end"] = aux.apply(lambda x: x[3]).astype(int)
            # df["new_start"] = df.apply(lambda row: row.start + row.start_clone, axis=1)
            # df["new_end"] = df.apply(lambda row: row.new_start + row.end_clone, axis=1)
            # df["lenght"] = df["new_end"] - df["new_start"]
            aux = aux.loc[aux["lenght"] > 30].reset_index(drop=True)
            aux = aux.loc[aux["lenght"] < 1200].reset_index(drop=True)
            aux = aux.sort_values(by=[0, "start"], ascending=[True, True])
            # aux[["ID", "chr", "new_start", "new_end", "lenght", "ave_depth", "start", "end", "geneid"]].to_csv(
            #     self.out + "_domainnot_merged.tab", header=None, sep="\t", index=False)
            aux[[0, "start", "end", "ID", "lenght", "ave_depth", 1, 2]].to_csv(
                self.out + "_domainnot_merged.tab", header=None, sep="\t", index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg105)
            sys.exit(1)
        else:
            self.filelog.write(msg106)
            return self.out + "_domainnot_merged.tab"

    def parsingoutput(self, outputdomainparsing, outputbedtoolscoverage):
        """
        parsing output file
        :param outputdomainparsing:
        :param outputbedtoolscoverage:
        :return:
        """
        try:
            df1 = pd.read_csv(outputdomainparsing, sep="\t", header=None,
                              names=["Chr", "CloneStart", "CloneEnd", "GeneID", "lenght", "ave_depth", "GeneStart",
                                     "GeneEnd"])
            df2 = pd.read_csv(outputbedtoolscoverage, sep="\t", header=None, index_col=False,
                              names=["ID", "start", "end", "GeneID", "dot1", "strand", "dot2", "description", "read_count"])
            df3 = pd.merge(df1, df2, on="GeneID")
            df3[["Chr", "CloneStart", "CloneEnd", "lenght", "GeneID", "strand", "read_count", "ave_depth",
                 "description"]].to_csv(
                self.out + "_definition.txt",
                header=True, sep="\t", index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg107)
            sys.exit(1)
        else:
            self.filelog.write(msg108)
            return self.out + "_definition.txt"

    def cleantempfile(self):
        """
        Remove temporany files.
        :return:
        """
        templistfile = ["_DepthCoverageBed.txt", "_BreadthCoverageBed.txt", ".bed", "_bedfixed.tab", "_readcov.bed",
                        "_maxdepth.txt",
                        "_percentilefiltered.bed", "_raw_domains.txt", "_domainnot_merged.tab"]
        for item in templistfile:
            if os.path.isfile(self.out + item):
                os.remove(self.out + item)

