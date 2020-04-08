from output_message import *
import optparse
import sys
import pandas as pd
import pybedtools
import traceback
import os

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
reference_opts.add_option('--log', action="store", dest="log", default=None,
                          help='Path/filelog.log.')
parser.add_option_group(reference_opts)

options, args = parser.parse_args()


class SubtractionProkaryotic(object):

    def __init__(self, optparseinstance):
        self.inputistance = optparseinstance
        #
        self.enrichedcontrol =  self.inputistance.enrichedcontrol
        self.enrichedselection = self.inputistance.enrichedselection
        self.overlap = self.inputistance.overlap
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
            if os.path.isfile(self.out + "_subtraction.log"):
                self.inputfilelog = self.out + "_subtraction.log"
            else:
                open(self.out + "_subtraction.log", 'a', 0)
                self.inputfilelog = self.out + "_subtraction.log"

    def filelogstdoutwrite(self, msg):
        """
        Write information about script esecution
        :param msg:
        :return:
        """
        self.inputfilelogopen = open(str(self.inputfilelog), 'a', 0)
        self.inputfilelogopen.write(msg)

    def filelogerrorwrite(self, msg):
        """
        Write error message
        :param msg:
        :return:
        """
        self.inputfilelogopen = open(str(self.inputfilelog), 'a', 0)
        self.inputfilelogopen.write(traceback.format_exc())
        self.inputfilelogopen.write(msg)
        sys.exit(1)

    def intersectresults(self, controlinput, selectioninput, overlap):
        """
        Function for reverse intersection between control and selection
        :param controlinput:
        :param selectioninput:
        :param overlap:
        :return:
        """
        try:
            control = pybedtools.BedTool(controlinput)
            selection = pybedtools.BedTool(selectioninput)
            res = selection.intersect(control, f=overlap, v=True, header=True)
            df1 = pd.read_table(res.fn, sep="\t", header=0)
            df1.to_csv(self.out + '_subtraction.txt', sep="\t", header=True, index=False)
        except traceback:
            self.filelogerrorwrite(msg153)
        else:
            self.filelogstdoutwrite(msg154)
            return self.out + '_subtraction.txt'


if __name__ == '__main__':
    DictSubtraction = dict()
    ClassSubtraction = SubtractionProkaryotic(optparseinstance=options)
    DictSubtraction["bedtoolsintersect"] = ClassSubtraction.intersectresults(
        ClassSubtraction.inputistance.enrichedcontrol,
        ClassSubtraction.enrichedselection,
        ClassSubtraction.overlap)
