from output_message import *
import optparse
import pandas as pd
import pybedtools
import os
import subprocess
import traceback
import sys

parser = optparse.OptionParser(usage='python %prog Enrichment Prokaryotic', version='1.0',)
parser.add_option('--outputcontrol', action="store", dest="outputcontrol", default=None,
                  help='TAB file derived from pyinteraseq_domains_definition.py'
                       ' with suffix _mapping.tab')
parser.add_option('--outputarget', action="store", dest="outputarget", default=None,
                  help='TAB file derived from pyinteraseq_domains_definition.py'
                       'with suffix _mapping.tab ',)

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
        self.bedfile = None
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

    def parsing_input(self, inputfile, prefixout):
        self.inputfilelogopen = open(str(self.inputfilelog), 'a')
        try:
            dfa = pd.read_csv(inputfile, index_col=False, header=0, sep='\t', names=[
                "ID", "domain_start", "domain_end", "length", "read_count", "strand", "average_depth", "chr",
                "geneID", "start", "end_x", "description"])
            dfa[["ID", "domain_start", "domain_end", "geneID", "read_count", "strand"]].to_csv(
                self.out + prefixout + '.bed', header=None,
                sep='\t', index=False)
        except traceback:
            self.filelogstdoutwrite(msg109)
        else:
            self.filelogstdoutwrite(msg110)
            return self.out + prefixout + '.bed'

    def commondomain(self, target, control):
        self.inputfilelogopen = open(str(self.inputfilelog), 'a')
        try:
            a = pybedtools.BedTool(target)
            b = pybedtools.BedTool(control)
            aintersectb = a.intersect(b, wo=True, r=True, f=0.5)
            aintersectb.moveto(self.out + "_intersect_commondomains.bed")
        except traceback:
            self.filelogstdoutwrite(msg111)
        else:
            self.filelogstdoutwrite(msg112)
            return self.out + "_intersect_commondomains.bed"

    def uniquedomain(self, target, control):
        self.inputfilelogopen = open(str(self.inputfilelog), 'a')
        try:
            a = pybedtools.BedTool(target)
            b = pybedtools.BedTool(control)
            aintersectb = a.intersect(b, v=True, r=True, f=0.5)
            aintersectb.moveto(self.out + "_intersect_uniquedomains.bed")
        except traceback:
            self.filelogstdoutwrite(msg113)
        else:
            self.filelogstdoutwrite(msg114)
            return self.out + "_intersect_uniquedomains.bed"

    def parserforedger(self, filecommondomains):
        self.inputfilelogopen = open(str(self.inputfilelog), 'a')
        try:
            dfa = pd.read_csv(filecommondomains, index_col=False, header=None, sep='\t',
                              names=["ID_t", "domain_start_t", "domain_end_t", "geneID_t", "read_count_t", "strand_t",
                                     "ID_c", "domain_start_c", "domain_end_c", "geneID_c", "read_count_c", "strand_c",
                                     "length"])
            dfa[["ID_t", "domain_start_t", "domain_end_t", "ID_t", "read_count_t", "strand_t"]].to_csv(
                self.out + '_target.bed', header=None, sep='\t', index=False)
            dfa[["ID_c", "domain_start_c", "domain_end_c", "ID_c", "read_count_c", "strand_c"]].to_csv(
                self.out + '_control.bed', header=None, sep='\t', index=False)
        except traceback:
            self.filelogstdoutwrite(msg115)
        else:
            self.filelogstdoutwrite(msg116)
            return self.out + '_target.bed', self.out + '_control.bed'

    def edger(self, parserforedger):
        self.inputfilelogopen = open(str(self.inputfilelog), 'a')
        fnull = open(os.devnull, 'w')
        try:
            subprocess.check_call(['Rscript', '--vanilla',
                                   self.path_edgeR, parserforedger[1], parserforedger[0], 'genomic', self.outputid,
                                   self.outputfolder], stdout=fnull, stderr=fnull)
        except subprocess.CalledProcessError:
            self.filelogstdoutwrite(msg117)
        except traceback:
            self.filelogstdoutwrite(msg117)
        else:
            self.filelogstdoutwrite(msg118)
            return self.out + ".sign_genes_adjpvalue_0.05.txt"

    def parsingoutputcommon(self, edgeroutput, inputtargetfile):
        self.inputfilelogopen = open(str(self.inputfilelog), 'a')
        try:
            dfa = pd.read_csv(edgeroutput, index_col=False, header=0, sep='\t',
                              names=['ID', 'chr', 'start', 'end', 'test_scripttarget_raw', 'test_scriptcondition_raw',
                                     'test_scripttarget_norm', 'test_scriptcondition_norm', 'logFC', 'logCPM', 'PValue',
                                     'AdjPValue'])
            dfb = dfa.loc[lambda df: df.logFC > 0, :]
            dfc = dfb[['ID', 'start', 'end', 'logFC', 'PValue', 'AdjPValue']]
            dfd = pd.read_csv(inputtargetfile, index_col=False, header=None, sep='\t',
                              names=["ID", "domain_start", "domain_end", "length", "read_count",
                                     "strand", "average_depth", "chr", "geneID", "start", "end_x", "description"])
            dff = pd.merge(dfc, dfd, on='ID')
            dff.rename(columns={'start_y': 'start', 'end_x': 'end'}, inplace=True)
            dff[['ID', 'domain_start', 'domain_end', "length", "chr", "geneID", "start", "end", "strand",
                 "logFC", "PValue", "AdjPValue", "description"]].to_csv(
                self.out + '_common_intervals.bed', header=True, sep='\t', index=False)
        except traceback:
            self.filelogstdoutwrite(msg119)
        else:
            self.filelogstdoutwrite(msg120)
            return self.out + '_common_intervals.bed'

    def parsingoutputunique(self, uniquedomains, inputcontrolfile):
        self.inputfilelogopen = open(str(self.inputfilelog), 'a')
        try:
            dfa = pd.read_csv(uniquedomains, index_col=False, header=None, sep='\t',
                              names=['ID', 'start', 'end', 'GeneID', 'read', 'strand'])
            dfb = pd.read_csv(inputcontrolfile, index_col=False, header=None, sep='\t',
                              names=["ID", "domain_start", "domain_end", "length", "read_count",
                                     "strand", "average_depth", "chr", "geneID", "start", "end_x", "description"])
            dfc = pd.merge(dfa, dfb, on='ID')
            dfc.rename(columns={'start_y': 'start', 'end_x': 'end', "strand_x": "strand"}, inplace=True)
            dfc[['ID', 'domain_start', 'domain_end', "length", "chr", "geneID", "start", "end", "strand", "description"]].to_csv(
                self.out + '_unique_intervals.bed', header=True, sep='\t', index=False)
        except traceback:
            self.filelogstdoutwrite(msg121)
        else:
            self.filelogstdoutwrite(msg122)
            return self.out + '_unique_intervals.bed'

    def cleantempfile(self):
        """
        Remove temporany files.
        :return:
        """
        templistfile = ["_intersect_commondomains.bed", "_intersect_uniquedomains.bed",
                        "_target.bed", "_control.bed", "_edger.txt", ".sign_genes_adjpvalue_0.05.txt", ""]
        for item in templistfile:
            if os.path.isfile(self.out + item):
                os.remove(self.out + item)


if __name__ == '__main__':
    DictEnrichment = dict()
    ClassEnrichment = EnrichmentProkaryotic(optparseinstance=options)
    DictEnrichment["controlparsed"] = ClassEnrichment.parsing_input(options.outputcontrol, "_control")
    DictEnrichment["targetparsed"] = ClassEnrichment.parsing_input(options.outputarget, "_target")
    DictEnrichment["common"] = ClassEnrichment.commondomain(DictEnrichment["targetparsed"],
                                                            DictEnrichment["controlparsed"])
    DictEnrichment["unique"] = ClassEnrichment.uniquedomain(DictEnrichment["targetparsed"],
                                                            DictEnrichment["controlparsed"])
    DictEnrichment["rfiles"] = ClassEnrichment.parserforedger(DictEnrichment["common"])
    DictEnrichment["edger"] = ClassEnrichment.edger(DictEnrichment["rfiles"])
    DictEnrichment["parsingcommon"] = ClassEnrichment.parsingoutputcommon(DictEnrichment["edger"],
                                                                          options.outputarget)
    DictEnrichment["parsingunique"] = ClassEnrichment.parsingoutputunique(DictEnrichment["unique"],
                                                                          options.outputcontrol)
    ClassEnrichment.cleantempfile()
