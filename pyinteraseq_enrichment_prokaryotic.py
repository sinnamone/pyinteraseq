from pyinteraseq_genomefileparsing import *
import optparse
import pandas as pd
import pybedtools
import os
import subprocess
from output_message import *

parser = optparse.OptionParser(usage='python %prog Enrichment Prokaryotic', version='1.0',)
parser.add_option('--blastnoutputgenomic', action="store", dest="blastnoutputgenomic", default=None,
                  help='Read dataset concatenate derived from pyinteraseq_mapping.py')
parser.add_option('--blastnoutputarget', action="store", dest="blastnoutputarget", default=None,
                  help='Read dataset concatenate derived from pyinteraseq_mapping.py',)
parser.add_option('--domainstarget', action="store", dest="domainstarget", default=None,
                  help='Output target derived from pyinteraseq_mapping.py',)

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
reference_opts.add_option('--thread', action="store", dest="thread", default='2',
                          help='Number of thread.')
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
            sys.stdout.write(msg9)
            sys.exit(0)
        #
        if self.inputistance.outputid is not None:
            self.outputid = self.inputistance.outputid
        else:
            sys.stdout.write(msg10)
            sys.exit(0)
        #
        self.out = self.outputfolder + self.outputid

    def blastninputparsing(self, fileinput, idexit):
        """
        :param fileinput: input file is created by pytinteraseq_domains_definition, and it has the following suffix
        \'_filtered_paired_complete.tab\'
        :param idexit: suffix to append to output file
        :return: bed 4 columns, chr, start, end, seqid related to mapping reads
        """
        df = pd.read_csv(fileinput, sep="\t", header=None)
        df1 = df[[1, 8, 9, 0]]
        df1 = df1.rename(columns={1: 0, 8: 1, 9: 2, 0: 3})
        df2 = df1[df1[1] > df1[2]].reset_index(drop=True)
        df3 = df1[df1[1] < df1[2]].reset_index(drop=True)
        df2 = df2[[0, 2, 1, 3]].reset_index(drop=True)
        df2 = df2.rename(columns={2: 1, 1: 2})
        df4 = df2.append(df3)
        df4.sort_values(by=[1]).to_csv(self.out + idexit + '.tab', sep="\t",
                                       header=None, index=False)
        return self.out + idexit + '.tab'

    def domainsdefinitionparsing(self, fileinput, idexit):
        """

        :return:
        """
        self.df1 = pd.read_csv(fileinput, sep="\t", header=None)
        self.df1[[0, 1, 2, 6, 8, 7]].sort_values(by=[1]).to_csv(self.out + idexit + '.tab',
                                                                sep="\t", header=None, index=False)
        return self.out + idexit + '.tab'

    def bedtoolscoveragebackgroung(self, outputdomainsparsed, outputblastnparsed, idex):
        """
        perferm bedtools coverage to get read count
        :param outputdomainsparsed:
        :param outputblastnparsed:
        :param idex:
        :return:
        """
        feata = pybedtools.BedTool(outputdomainsparsed)
        featb = pybedtools.BedTool(outputblastnparsed)
        c = feata.coverage(featb, counts=True)
        df1 = pd.read_table(c.fn, names=['chr', 'clonestart', 'cloneend',
                                         'chr2', 'start', 'end', 'geneid',
                                         'cog', 'strand', 'genename', 'description', 'clonelength'])
        df1.to_csv(self.out + idex + '.bed', sep="\t", header=None, index=False)
        return self.out + idex + '.bed'

    def edgerinputformatparsing(self, outputbedtoolscoverage, idex):
        """

        :param outputbedtoolscoverage:
        :param idex:
        :return:
        """
        self.df1 = pd.read_csv(outputbedtoolscoverage, sep="\t", header=None)
        # add column with dataframe index
        self.df1['index1'] = self.df1.index
        # merge gene id with index in order to create unique geneid
        self.df1['geneid'] = self.df1[3] + '_' + self.df1["index1"].astype(str)
        # write output
        self.df1[[0, 1, 2, 'geneid', 6, 5]].to_csv(self.out + idex + '.bed',
                                                   sep="\t", header=None, index=False)
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
        subprocess.check_call(['/Library/Frameworks/R.framework/Versions/3.2/Resources/bin/Rscript', '--vanilla',
                               self.path_edgeR, backgroundedgeparsed, targetedgeparsed, 'genomic', idex, outpath])
        return outpath + 'genomic_' + idex + '.sign_genes_adjpvalue_0.05.txt'

    def edgeroutparser(self, edgeroutput, outputdomaindetection):
        """

        :param edgeroutput:
        :param outputdomaindetection:
        :return:
        """
        df1 = pd.read_csv(edgeroutput, sep="\t", header=0)
        df2 = df1.loc[(df1['logFC'] > 0.0)]
        df3 = pd.read_csv(outputdomaindetection, sep="\t", header=None,
                          names=['chr', 'clonestart', 'cloneend', 'clonelength', 'start',
                                 'end', 'geneid','strand', 'genename', 'description', 'nseq'])
        df4 = pd.merge(df3, df2, right_on='start', left_on='clonestart')
        df4[['chr_x', 'clonestart', 'cloneend', 'clonelength', 'start_x', 'end_x', 'geneid', 'logFC', 'PValue',
             'AdjPValue', 'strand', 'genename', 'description', 'nseq']].to_csv(
            self.out+'_enrichment.txt', sep="\t", header=None, index=False)
        return self.out+'_enrichment.txt'


if __name__ == '__main__':
    DictEnrichment = dict()
    # parsing input file, blastn output filtered both for background and genomic
    DictEnrichment["blastnbackgroundparsed"] = EnrichmentProkaryotic(optparseinstance=options).blastninputparsing(
        fileinput=options.blastnoutputgenomic, idexit='_blastnbackgroundparsed')
    DictEnrichment["blastntargetparsed"] = EnrichmentProkaryotic(optparseinstance=options).blastninputparsing(
        fileinput=options.blastnoutputarget, idexit='_blastntargetparsed')
    DictEnrichment["outputtargetparsed"] = EnrichmentProkaryotic(optparseinstance=options).domainsdefinitionparsing(
        fileinput=options.domainstarget, idexit='_outputtargetparsed')
    DictEnrichment["bedcoveragebackground"] = EnrichmentProkaryotic(optparseinstance=options).bedtoolscoveragebackgroung(
        DictEnrichment["outputtargetparsed"], DictEnrichment["blastnbackgroundparsed"], '_back')
    DictEnrichment["bedcoveragetarget"] = EnrichmentProkaryotic(optparseinstance=options).bedtoolscoveragebackgroung(
        DictEnrichment["outputtargetparsed"], DictEnrichment["blastntargetparsed"], '_target')
    DictEnrichment["bedcoveragebackgroundready"] = EnrichmentProkaryotic(optparseinstance=options).\
        edgerinputformatparsing(DictEnrichment["bedcoveragebackground"], '_backedgeready')
    DictEnrichment["bedcoveragetargetready"] = EnrichmentProkaryotic(optparseinstance=options).edgerinputformatparsing(
        DictEnrichment["bedcoveragetarget"], '_targetedgeready')
    DictEnrichment["edgertable"] = EnrichmentProkaryotic(optparseinstance=options).edger(
        DictEnrichment["bedcoveragebackgroundready"], DictEnrichment["bedcoveragetargetready"],
        options.outputid, options.outputfolder)
    EnrichmentProkaryotic(optparseinstance=options).edgeroutparser(
        DictEnrichment["edgertable"], options.domainstarget)

