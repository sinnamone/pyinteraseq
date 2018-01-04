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


#
#     def cleantempfile(self):
#         """
#         Remove temporany files.
#         :return:
#         """
#         db1 = str(self.out + ".nsq")
#         db2 = str(self.out + ".nin")
#         db3 = str(self.out + ".nhr")
#         os.remove(db1)
#         os.remove(db2)
#         os.remove(db3)
#         templistfile = ["_clones.tab","_newid.tab","_renamed.fasta","_otus_most_abundant.fa","_clean.fasta","_blastnclones.tab","_blastclonesparsed.bed",
#                         "_bedfixed.tab","_clonesannotated.bed","_clonesannotatedfiltered.bed","_cluster_count.txt","_blastnclonescounted.bed","_blastnclonescountedfiltered.bed",
#                         "_blastnclonesmerge.bed","fasta_seq.tab","_clean.tab","_clonestabular.tab","_clonesdescription.bed"]
#         for item in templistfile:
#             if os.path.isfile(self.out + item):
#                 os.remove(self.out + item)
#         pickedlist = ["_renamed_otus.log","_renamed_otus.txt","_renamed_clusters.uc"]
#         for item in pickedlist:
#             if os.path.isfile(self.out + '_picked/' + self.outputid + item):
#                 os.remove(self.out + '_picked/' + self.outputid + item)
#         os.rmdir(self.out + '_picked/')

class DomainsDefinition(InputCheckDomainDefinition):

    def __init__(self, optparseinstance):
        InputCheckDomainDefinition.__init__(self, optparseinstance)

    def depthcoverage(self):
        self.filelog = self.logopen()
        try:
            with open(self.out + '_DepthCoverageBed.txt', 'w') as b:
                subprocess.check_call(['bedtools', 'genomecov','-dz', '-ibam', self.bamfile, '-g', self.genome],stdout=b)
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
        self.filelog = self.logopen()
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
        :return:
        """
        self.filelog = self.logopen()
        try:
            bedfile = pybedtools.BedTool(self.bamfile).bam_to_bed()
            bedfile = pybedtools.BedTool(self.bamfile).bam_to_bed().saveas(self.out + '.bed')
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
        self.filelog = self.logopen()
        try:
            self.df1 = pd.read_csv(convertedbam, sep="\t", header=None,
                                   names=["ID", "start_clone", "end_clone", "id_clone", "score", "strand"])
            self.aux = self.df1["ID"].apply(lambda x: x.split('_'))
            self.df1["start"] = self.aux.apply(lambda x: x[2]).astype(int)
            self.df1["end"] = self.aux.apply(lambda x: x[3]).astype(int)
            self.df1["new_start"] = self.df1.apply(lambda row: row.start + row.start_clone, axis=1)
            self.df1["new_end"] = self.df1.apply(lambda row: row.new_start + row.end_clone, axis=1)
            self.df1[["ID", "new_start", "new_end", "id_clone", "score", "strand"]].to_csv(self.out + "_bedfixed.tab", header=None, sep="\t", index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg139)
            sys.exit(1)
        else:
            self.filelog.write(msg140)
            return self.out + "_bedfixed.tab"

    def bedtoolscoverage(self, outputcoordinatesfix):
        try:
            bed = pybedtools.BedTool(outputcoordinatesfix)
            annotation = pybedtools.BedTool(self.annotation)
            cov = annotation.coverage(bed,counts=True).saveas(self.out + '_readcov.bed')
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg98)
            sys.exit(1)
        else:
            self.filelog.write(msg97)
            return self.out + '_readcov.bed'

    def groupbydepth(self, outputbreadthcoverage):
        try:
            df = pd.read_csv(outputbreadthcoverage, header=None, sep='\t', names=["chr","depth","numbasdepth","lenght","perccov"])
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
        try:
            df = pd.read_csv(outputstartenddefinition,sep="\t", header=None,
                             names=["ID", "start_clone", "end_clone", "ave_depth"])
            aux = df["ID"].apply(lambda x: x.split('_'))
            df["geneid"] = aux.apply(lambda x: x[0])
            df["chr"] = aux.apply(lambda x: x[1])
            df["start"] = aux.apply(lambda x: x[2]).astype(int)
            df["end"] = aux.apply(lambda x: x[3]).astype(int)
            df["new_start"] = df.apply(lambda row: row.start + row.start_clone, axis=1)
            df["new_end"] = df.apply(lambda row: row.new_start + row.end_clone, axis=1)
            df["lenght"] = df["new_end"] - df["new_start"]
            df = df.loc[df["lenght"] > 30].reset_index(drop=True)
            df = df.loc[df["lenght"] < 1200].reset_index(drop=True)
            df = df.sort_values(by=["chr","start"],ascending=[True,True])
            df[["ID","chr","new_start", "new_end","lenght", "ave_depth", "start", "end","geneid"]].to_csv(self.out + "_domainnot_merged.tab",
                                                                                   header=None, sep="\t", index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg105)
            sys.exit(1)
        else:
            self.filelog.write(msg106)
            return self.out + "_domainnot_merged.tab"

    def parsingoutput(self,outputdomainparsing,outputbedtoolscoverage):
        try:
            df1 = pd.read_csv(outputdomainparsing,sep="\t", header=None,
                              names=["ID","chr","domain_start", "domain_end","length", "average_depth", "start", "end","geneID"])
            df2 = pd.read_csv(outputbedtoolscoverage,sep="\t", header=None,
                              names=["ID","start", "end","geneID","score","strand","description","read_count"])
            df3 = pd.merge(df1,df2,on="ID")
            df3[["ID","domain_start", "domain_end","length","read_count","strand","average_depth","chr","geneID_x","start_x", "end_x","description"]].to_csv(
                self.out + "_domain_definition.tab",
                header=True, sep="\t", index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg107)
            sys.exit(1)
        else:
            self.filelog.write(msg108)
            return self.out + "_domain_definition.tab"


    # def clean(self, dest):
    #     try:
    #         suffixdelete = ["temp_1.bed", "temp_2.bed", "temp_3.bed", "temp_4.bed", "temp_5.bed", "temp_6.bed",
    #                         "temp_7.bed", "temp_8.bed", "temp_9.bed", "temp_10.bed"]
    #         for file in os.listdir(dest):
    #             for i in range(len(suffixdelete)):
    #                 if file.endswith(suffixdelete[i]):
    #                     os.remove(os.path.join(dest, file))
    #     except ValueError:
    #         sys.stdout.write('Error. Removing temporary files. Exit\n')
    #         sys.exit(0)
    #     else:
    #         sys.stdout.write('Removing temporary files. Complete.\n')
    #


