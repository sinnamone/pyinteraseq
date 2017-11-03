from pyinteraseq_genomefileparsing import *
import optparse
import pandas as pd
import pybedtools
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from rpy2.robjects import *
from output_message import *

parser = optparse.OptionParser(usage='python %prog Enrichment Prokariotyc', version='1.0',)
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

    # def edger(self, edgerinputformatparsingtarget, edgerinputformatparsingbackground):
    #     edgeR = importr('edgeR')
    #     r.assign('infile_cond1', '/Users/simone/output_test/SelA_versus_26695_targetedgeready.bed')
    #     r.assign('infile_cond2', '/Users/simone/output_test/SelA_versus_26695_backedgeready.bed')
    #     r.assign('name_cond1', 'target')
    #     r.assign('name_cond2', 'background')
    #     r.assign('output_destination', self.outputfolder)
    #     r.assign('threshold_fdr', 0.05)
    #     r('bed_cond1 =  read.table(infile_cond1,header=FALSE)')
    #     r('counts_cond1 = bed_cond1[,4:5]')
    #     r('colnames(counts_cond1) = c("ID",name_cond1)')
    #     r('counts_cond1 = counts_cond1[order(counts_cond1[,1]),]')
    #     r('bed_cond2 =  read.table(infile_cond2,header=FALSE)')
    #     r('counts_cond2 = bed_cond2[,4:5]')
    #     r('colnames(counts_cond2) = c("ID",name_cond2)')
    #     r('counts_cond2 = counts_cond2[order(counts_cond2[,1]),]')
    #     r('counts = merge(counts_cond1, counts_cond2, by="ID")')
    #     r('rownames(counts) = counts[, 1]')
    #     r('counts = counts[, c(2:3)]')
    #     r('group = factor(c(1, 2))')
    #     r('dge = DGEList(counts=counts, group=group)')
    #     res = r('y = calcNormFactors(dge)')
    #     print res
        # robjects.r('et2 = exactTest(y, dispersion="auto")')
        # res15 = robjects.r('counts = counts[order(rownames(counts)),]')
        # res16 = robjects.r('diff = et2$table')
        # res17 = robjects.r('diff = diff[order(rownames(diff)),]')
        # res18 = robjects.r('norm = counts')
        # res19 = robjects.r('norm[, 1] = counts[, 1]*(y$samples[name_cond1, "norm.factors"])')
        # res20 = robjects.r('norm[, 2] = counts[, 2]*(y$samples[name_cond2, "norm.factors"])')
        # res21 = robjects.r('output = cbind(counts, norm, diff)')
        # res22 = robjects.r('output = output[order(rownames(output)),]')
        # res23 = robjects.r('bed_cond1 = bed_cond1[order(bed_cond1[, 4]), ]')
        # res24 = robjects.r('bed_cond2 = bed_cond2[order(bed_cond2[, 4]), ]')
        # res25 = robjects.r('chr = bed_cond1[, 1]')
        # res26 = robjects.r('start = bed_cond1[, 2]')
        # res27 = robjects.r('end = bed_cond1[, 3]')
        # res28 = robjects.r('output = cbind(chr, start, end, output)')
        # res29 = robjects.r('output = output[order(output[, 1], output[, 2]), ]')
        # res30 = robjects.r('colnames(output)[4:5] = paste(colnames(output)[4:5], "_raw", sep="")')
        # res31 = robjects.r('colnames(output)[6:7] = paste(colnames(output)[6:7], "_norm", sep="")')
        # res32 = robjects.r('ID = rownames(output)')
        # res33 = robjects.r('AdjPValue = p.adjust(output$PValue, method = "fdr")')
        # res34 = robjects.r('output = cbind(ID, output, AdjPValue)')
        # res35 = robjects.r('write.table(output, file=paste(output_destination, name_cond1, name_cond1, "_", name_cond2, ".txt", sep=""),'
        #                    'quote=FALSE, sep="\t", row.names = FALSE)')
        #
        # res36 = robjects.r('output_sign = output[output[, "AdjPValue"] < threshold_fdr,]')
        # res37 = robjects.r('write.table(output_sign,file=paste(output_destination, name_cond1, "_", '
        #                    'name_cond2, ".sign_genes_adjpvalue_", threshold_fdr,".txt", sep=""), '
        #                    'quote=FALSE, sep="\t", row.names = FALSE)')

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
    # EnrichmentProkaryotic(optparseinstance=options).edger(
    #     DictEnrichment["bedcoveragetargetready"], DictEnrichment["bedcoveragebackgroundready"])
