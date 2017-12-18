import optparse
from output_message import *
import pybedtools
import sys
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib_venn import venn3
import traceback
import matplotlib
plt.switch_backend('agg')

parser = optparse.OptionParser(usage='python %prog Intersection Prokaryotic', version='1.0')
parser.add_option('--selection', action="append", dest='selections', default=[])
parser.add_option('--label', action="append", dest='labels', default=[])
parser.add_option('--outputfolder', action="store", dest="outputfolder",
                  default=None, help='Output folder.')
parser.add_option('--outputid', action="store", dest="outputid",
                  default=None, help='Output ID.')
parser.add_option('--log', action="store", dest="log", default=None, help='Path/filelog.log.')
options, args = parser.parse_args()

filelog = None
count = 0
interection = options.outputid + '_intersection.txt'
unique = options.outputid + '_unique.txt'
for i in range(len(options.selections)):
    count += 1
if options.log is None:
    filelog = options.outputfolder + options.outputid + '_intersection.log'
elif options.log is not None:
    filelog = options.log

if count == 0:
    filelogopen = open(str(filelog), 'a')
    filelogopen.write(msg85)
    sys.exit(1)
elif count == 1:
    filelogopen = open(str(filelog), 'a')
    filelogopen.write(msg86)
    sys.exit(1)
elif count == 2:
    try:
        filelogopen = open(str(filelog), 'a')
        a = pybedtools.BedTool(options.selections[0])
        b = pybedtools.BedTool(options.selections[1])
        # intersection between two list
        aintersectb = a.intersect(b, u=True)
        # unique for list A
        aunique = a.intersect(b, v=True)
        # unique for list B
        bunique = b.intersect(a, v=True)
        v = venn2(subsets=((a - b).count(), (a + b).count(), (b - a).count()),
                  set_labels=(options.labels[0], options.labels[1]))
        plt.savefig(options.outputfolder + options.outputid + '.pdf', format='pdf')
        # write files
        (a - b).moveto(options.outputfolder + options.outputid + '_'.join([options.labels[0], unique]))
        (b - a).moveto(options.outputfolder + options.outputid + '_'.join([options.labels[1], unique]))
        (a + b).moveto(options.outputfolder + options.outputid + '_'.join([options.labels[0], options.labels[1],
                                                                           interection]))
    except traceback:
        filelogopen = open(str(filelog), 'a')
        filelogopen.write(msg155)
        sys.exit(1)
    else:
        filelogopen.write(msg156)
        sys.exit(1)

elif count == 3:
    try:
        filelogopen = open(str(filelog), 'a')
        a = pybedtools.BedTool(options.selections[0])
        b = pybedtools.BedTool(options.selections[1])
        c = pybedtools.BedTool(options.selections[2])
        v = venn3(subsets=(
            (a - b - c).count(), (b - a - c).count(), (a + b - c).count(), (c - a - b).count(), (c + a - b).count(),
            (c + b - a).count(), (a + b + c).count()), set_labels=(options.labels[0], options.labels[1],
                                                                   options.labels[2]))
        plt.savefig(options.outputfolder + options.outputid + '.pdf', format='pdf')
        #
        (a + b).moveto(options.outputfolder
                       + '_'.join([options.labels[0], options.labels[1], interection]))
        (a + c).moveto(options.outputfolder
                       + '_'.join([options.labels[0], options.labels[2], interection]))
        (b + c).moveto(options.outputfolder
                       + '_'.join([options.labels[1], options.labels[2], interection]))
        (a + b + c).moveto(options.outputfolder
                           + '_'.join([options.labels[0], options.labels[1], options.labels[2] + '.txt']))
        (a - b - c).moveto(options.outputfolder
                           + '_'.join([options.labels[0], unique]))
        (b - a - c).moveto(options.outputfolder
                           + '_'.join([options.labels[1], unique]))
        (c - a - b).moveto(options.outputfolder
                           + '_'.join([options.labels[2], unique]))

    except traceback:
        filelogopen = open(str(filelog), 'a')
        filelogopen.write(msg155)
        sys.exit(1)
    else:
        filelogopen.write(msg156)
        sys.exit(1)

