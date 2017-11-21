from output_message import *
import optparse
import sys
import pandas as pd
import pybedtools

parser = optparse.OptionParser(usage='python %prog Subtraction Prokaryotic', version='1.0',)
parser.add_option('--enrichedcontrol', action="store", dest="enrichedcontrol", default=None,
                  help='')
parser.add_option('--enrichedselection', action="store", dest="enrichedselection", default=None,
                  help='')

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
reference_opts.add_option('--overlap', action="store", dest="overlap", default='0.50',
                          help='Minimum overlap required.')
reference_opts.add_option('--thread', action="store", dest="thread", default='2',
                          help='Number of thread.')
parser.add_option_group(reference_opts)

options, args = parser.parse_args()


class SubtractionProkaryotic(object):

    def __init__(self, optparseinstance):
        self.inputistance = optparseinstance
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

    def intersectresults(self, controlinput, selectioninput, overlap):
        """
        Function for reverse intersection between control and selection
        :param controlinput:
        :param selectioninput:
        :param overlap:
        :return:
        """
        control = pybedtools.BedTool(controlinput)
        selection = pybedtools.BedTool(selectioninput)
        res = selection.intersect(control, f=overlap, v=True)
        df1 = pd.read_table(res.fn, sep="\t", header=None)
        df1.to_csv(self.out + '_subtraction.txt', sep="\t", header=None, index=False)
        return self.out + '_subtraction.txt'


if __name__ == '__main__':
    DictSubtraction = dict()
    DictSubtraction["bedtoolsintersect"] = SubtractionProkaryotic(optparseinstance=options).intersectresults(
        options.enrichedcontrol, options.enrichedselection, options.overlap)
