from pyinteraseq_genomefileparsing import *
from output_message import *
import optparse
import pandas as pd
import pybedtools
import os
import subprocess
import traceback

parser = optparse.OptionParser(usage='python %prog Enrichment Prokaryotic', version='1.0',)
parser.add_option('--blastnoutputgenomic', action="store", dest="blastnoutputgenomic", default=None,
                  help='Read dataset concatenate derived from pyinteraseq_mapping.py'
                       ' with suffix _mapping.tab')
parser.add_option('--blastnoutputarget', action="store", dest="blastnoutputarget", default=None,
                  help='Read dataset concatenate derived from pyinteraseq_mapping.py '
                       'with suffix _mapping.tab ',)
parser.add_option('--domainstarget', action="store", dest="domainstarget", default=None,
                  help='Output target derived from pyinteraseq_domain_definition.py'
                       ' with suffix _domain_detection.tab ',)

query_opts = optparse.OptionGroup(
    parser, 'Output Options',
    'Options for the output destionation and name.',
    )
query_opts.add_option('--outputfolder', action="store", dest="outputfolder", default=None,
                      help='Output folder.')
query_opts.add_option('--outputid', action="store", dest="outputid", default=None,
                      help='Output ID.')

parser.add_option_group(query_opts)


reference_opts = optparse.OptionGroup(
    parser, 'Advanced Options',
    'Options for advanced analysis.',
    )
reference_opts.add_option('--log', action="store", dest="log", default=None,
                          help='Path/filelog.log.')
parser.add_option_group(reference_opts)

options, args = parser.parse_args()


class EnrichmentProkaryotic(object):

    def __init__(self, optparseinstance):
        self.inputistance = optparseinstance
        self.df1 = None
        self.path_edgeR = os.path.dirname(os.path.realpath(__file__)) + '/diff_edgeR.Rscript'
        #
        if self.inputistance.outputfolder is not None:
            if self.inputistance.outputfolder.endswith('/') is True:
                self.outputfolder = self.inputistance.outputfolder
            else:
                self.outputfolder = self.inputistance.outputfolder + '/'
        else:
            self.filelogerrorwrite(msg9)
        #
        if self.inputistance.outputid is not None:
            self.outputid = self.inputistance.outputid
        else:
            self.filelogerrorwrite(msg10)
        #
        self.out = self.outputfolder + self.outputid
        self.inputfilelogopen = None
        self.inputfilelog = self.inputistance.log
        if self.inputfilelog is None:
            if os.path.isfile(self.out + "_enrichment.log"):
                self.inputfilelog = self.out + "_enrichment.log"
            else:
                open(self.out + "_enrichment.log", 'a')
                self.inputfilelog = self.out + "_enrichment.log"

    def filelogstdoutwrite(self, msg):
        """
        Write information about script esecution
        :param msg:
        :return:
        """
        self.inputfilelogopen = open(str(self.inputfilelog), 'a')
        self.inputfilelogopen.write(msg)

    def filelogerrorwrite(self, msg):
        """
        Write error message
        :param msg:
        :return:
        """

        self.inputfilelogopen = open(str(self.inputfilelog), 'a')
        self.inputfilelogopen.write(traceback.format_exc())
        self.inputfilelogopen.write(msg)
        sys.exit(1)

    def blastninputparsing(self, fileinput, idexit):
        """
        :param fileinput: input file is created by pytinteraseq_mapping, and it has the following suffix
        \'_output_mapping_step.tab\'
        :param idexit: suffix to append to output file
        :return: bed 4 columns, chr, start, end, seqid related to mapping reads
        """
        try:
            df = pd.read_csv(fileinput, sep="\t", header=0)
            df.columns = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
            df1 = df[[1, 8, 9, 0]]
            df1 = df1.rename(columns={1: 0, 8: 1, 9: 2, 0: 3})
            df2 = df1[df1[1] > df1[2]].reset_index(drop=True)
            df3 = df1[df1[1] < df1[2]].reset_index(drop=True)
            df2 = df2[[0, 2, 1, 3]].reset_index(drop=True)
            df2 = df2.rename(columns={2: 1, 1: 2})
            df4 = df2.append(df3)
            df4.sort_values(by=[1]).to_csv(self.out + idexit + '.tab', sep="\t",
                                           header=None, index=False)
        except traceback:
            self.filelogerrorwrite(msg139)
        else:
            self.filelogstdoutwrite(msg140)
            return self.out + idexit + '.tab'

    def domainsdefinitionparsing(self, fileinput, idexit):
        """
        Parser function of pyinteraseq_domain_definition.py output
        :return:
        """
        try:
            self.df1 = pd.read_csv(fileinput, sep="\t", header=0)
            self.df1.columns = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            self.df1[[0, 1, 2, 6, 8, 7]].sort_values(by=[1]).to_csv(self.out + idexit + '.tab',
                                                                    sep="\t", header=None, index=False)
        except traceback:
            self.filelogerrorwrite(msg141)
        else:
            self.filelogstdoutwrite(msg142)
            return self.out + idexit + '.tab'

    def bedtoolscoveragebackgroung(self, outputdomainsparsed, outputblastnparsed, idex):
        """
        perferm bedtools coverage to get read count
        :param outputdomainsparsed:
        :param outputblastnparsed:
        :param idex:
        :return:
        """
        try:
            feata = pybedtools.BedTool(outputdomainsparsed)
            featb = pybedtools.BedTool(outputblastnparsed)
            c = feata.coverage(featb, counts=True)
            df1 = pd.read_table(c.fn, names=['chr', 'clonestart', 'cloneend',
                                             'chr2', 'start', 'end', 'geneid',
                                             'cog', 'strand', 'genename', 'description', 'clonelength'])
            df1.to_csv(self.out + idex + '.bed', sep="\t", header=None, index=False)
        except traceback:
            self.filelogerrorwrite(msg143)
        else:
            self.filelogstdoutwrite(msg144)
            return self.out + idex + '.bed'

    def edgerinputformatparsing(self, outputbedtoolscoverage, idex):
        """
        Parser output bedtools coverage
        :param outputbedtoolscoverage:
        :param idex:
        :return: Bed6
        """
        try:
            self.df1 = pd.read_csv(outputbedtoolscoverage, sep="\t", header=None)
            # filter
            self.df1 = self.df1.loc[self.df1[3] != "."]
            # add column with dataframe index
            self.df1['index1'] = self.df1.index
            # merge gene id with index in order to create unique geneid
            self.df1['geneid'] = self.df1[3] + '_' + self.df1["index1"].astype(str)
            # write output
            self.df1[[0, 1, 2, 'geneid', 6, 5]].to_csv(self.out + idex + '.bed',
                                                       sep="\t", header=None, index=False)
        except traceback:
            self.filelogerrorwrite(msg147)
        else:
            self.filelogstdoutwrite(msg148)
            return self.out + idex + '.bed'

    def edger(self, backgroundedgeparsed, targetedgeparsed, idex, outpath):
        """
        edgeR function
        :param backgroundedgeparsed:
        :param targetedgeparsed:
        :param idex:
        :param outpath:
        :return:
        """
        self.inputfilelogopen = open(str(self.inputfilelog), 'a')
        fnull = open(os.devnull, 'w')
        try:
            subprocess.check_call(['Rscript', '--vanilla',
                                   self.path_edgeR, backgroundedgeparsed, targetedgeparsed, 'genomic', idex, outpath],
                                  stderr=fnull, stdout=fnull)
        except subprocess.CalledProcessError:
            self.inputfilelog.write(msg149)
        except traceback:
            self.filelogerrorwrite(msg149)
        else:
            self.filelogstdoutwrite(msg150)
            return outpath + 'genomic_' + idex + '.sign_genes_adjpvalue_0.05.txt'

    def edgeroutparser(self, edgeroutput, outputdomaindetection):
        """
        Parser and filtering function of edgeR test
        :param edgeroutput:
        :param outputdomaindetection:
        :return:
        """
        try:
            df1 = pd.read_csv(edgeroutput, sep="\t", header=0)
            df2 = df1.loc[(df1['logFC'] > 0.0)]
            df3 = pd.read_csv(outputdomaindetection, sep="\t", header=0,
                              names=['chr', 'clonestart', 'cloneend', 'clonelength', 'start',
                                     'end', 'geneid', 'strand', 'genename', 'description', 'nseq'])
            df4 = pd.merge(df3, df2, right_on='start', left_on='clonestart')
            df5 = df4[['chr_x', 'clonestart', 'cloneend', 'clonelength', 'start_x', 'end_x', 'geneid', 'logFC', 'PValue',
                       'AdjPValue', 'strand', 'genename', 'description', 'nseq']]
            df5.columns = ["Chr", "CloneStart", "CloneEnd", "CloneLength", "Start", "End", "GeneID", "logFC",
                           "PValue", "AdjPValue", "Strand", "GeneName", "Description", "NuclSeq"]
            df5.to_csv(self.out + '_enrichment.txt', sep="\t", header=True, index=False)
        except traceback:
            self.filelogerrorwrite(msg151)
        else:
            self.filelogstdoutwrite(msg152)
            return self.out+'_enrichment.txt'

    def cleantempfiles(self):
        """
        clean temporany files
        :return:
        """
        templistfile = ["_blastnbackgroundparsed.tab", "_blastntargetparsed.tab", "_outputtargetparsed.tab",
                        "_back.bed", "_target.bed",
                        "_backedgeready.bed", "_targetedgeready.bed"]
        for item in templistfile:
            if os.path.isfile(self.out + item):
                os.remove(self.out + item)
        if os.path.isfile(self.outputfolder + "genomicgenomic_" + self.outputid + ".txt"):
            os.remove(self.outputfolder + "genomicgenomic_" + self.outputid + ".txt")
        if os.path.isfile(self.outputfolder + "genomic_" + self.outputid + ".sign_genes_adjpvalue_0.05.txt"):
            os.remove(self.outputfolder + "genomic_" + self.outputid + ".sign_genes_adjpvalue_0.05.txt")


if __name__ == '__main__':
    DictEnrichment = dict()
    ClassEnrichment = EnrichmentProkaryotic(optparseinstance=options)
    # parsing input file, blastn output filtered both for background and genomic
    DictEnrichment["blastnbackgroundparsed"] = ClassEnrichment.blastninputparsing(
        fileinput=options.blastnoutputgenomic, idexit='_blastnbackgroundparsed')
    DictEnrichment["blastntargetparsed"] = ClassEnrichment.blastninputparsing(
        fileinput=options.blastnoutputarget, idexit='_blastntargetparsed')
    DictEnrichment["outputtargetparsed"] = ClassEnrichment.domainsdefinitionparsing(
        fileinput=options.domainstarget, idexit='_outputtargetparsed')
    DictEnrichment["bedcoveragebackground"] = ClassEnrichment.bedtoolscoveragebackgroung(
        DictEnrichment["outputtargetparsed"], DictEnrichment["blastnbackgroundparsed"], '_back')
    DictEnrichment["bedcoveragetarget"] = ClassEnrichment.bedtoolscoveragebackgroung(
        DictEnrichment["outputtargetparsed"], DictEnrichment["blastntargetparsed"], '_target')
    DictEnrichment["bedcoveragebackgroundready"] = ClassEnrichment.edgerinputformatparsing(
        DictEnrichment["bedcoveragebackground"], '_backedgeready')
    DictEnrichment["bedcoveragetargetready"] = ClassEnrichment.edgerinputformatparsing(
        DictEnrichment["bedcoveragetarget"], '_targetedgeready')
    DictEnrichment["edgertable"] = ClassEnrichment.edger(
        DictEnrichment["bedcoveragebackgroundready"], DictEnrichment["bedcoveragetargetready"],
        options.outputid, options.outputfolder)
    DictEnrichment["outputtab"] = ClassEnrichment.edgeroutparser(DictEnrichment["edgertable"], options.domainstarget)
    ClassEnrichment.cleantempfiles()
