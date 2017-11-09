import optparse
from output_message import *
import pybedtools
import sys
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib_venn import venn3

parser = optparse.OptionParser(usage='python %prog Intersection Prokaryotic', version='1.0')
parser.add_option('-i', action="append", dest='inputs', default=[])
parser.add_option('-l', action="append", dest='labels', default=[])
parser.add_option('--outputfolder', action="store", dest="outputfolder",
                  default=None, help='Output folder.')
parser.add_option('--outputid', action="store", dest="outputid",
                  default=None, help='Output ID.')
options, args = parser.parse_args()


count = 0
for i in range(len(options.inputs)):
    count += 1


if count == 0:
    sys.stdout.write(msg85)
    sys.exit(0)
elif count == 1:
    sys.stdout.write(msg86)
    sys.exit(0)
elif count == 2:
    a = pybedtools.BedTool(options.inputs[0])
    b = pybedtools.BedTool(options.inputs[1])
    # intersection between two list
    aintersectb = a.intersect(b, u=True)
    # unique for list A
    aunique = a.intersect(b, v=True)
    # unique for list B
    bunique = b.intersect(a, v=True)
    v = venn2(subsets=((a-b).count(), (a+b).count(), (b-a).count()), set_labels=(options.labels[0], options.labels[1]))
    plt.savefig(options.outputfolder + options.outputid + '.pdf', format='pdf')
    # write files
    (a-b).moveto(options.outputfolder + options.outputid + '_'.join([options.labels[0], 'unique.bed']))
    (b-a).moveto(options.outputfolder + options.outputid + '_'.join([options.labels[1], 'unique.bed']))
    (a+b).moveto(options.outputfolder + options.outputid + '_'.join([options.labels[0], options.labels[1],
                                                                     'intersection.bed']))
elif count == 3:
    a = pybedtools.BedTool(options.inputs[0])
    b = pybedtools.BedTool(options.inputs[1])
    c = pybedtools.BedTool(options.inputs[2])
    v = venn3(subsets=((a-b-c).count(), (b-a-c).count(), (a+b-c).count(), (c-a-b).count(), (c+a-b).count(),
                       (c+b-a).count(), (a+b+c).count()), set_labels=(options.labels[0], options.labels[1],
                                                                      options.labels[2]))
    plt.savefig(options.outputfolder + options.outputid + '.pdf', format='pdf')
    #
    (a+b).moveto(options.outputfolder + options.outputid
                 + '_'.join([options.labels[0], options.labels[1], 'intersection.bed']))
    (a+c).moveto(options.outputfolder + options.outputid
                 + '_'.join([options.labels[0], options.labels[2], 'intersection.bed']))
    (b+c).moveto(options.outputfolder + options.outputid
                 + '_'.join([options.labels[1], options.labels[2], 'intersection.bed']))
    (a+b+c).moveto(options.outputfolder + options.outputid
                   + '_'.join([options.labels[0], options.labels[1], options.labels[2] + '.bed']))
    (a-b-c).moveto(options.outputfolder + options.outputid
                    + '_'.join([options.labels[0], 'unique.bed']))
    (b-a-c).moveto(options.outputfolder + options.outputid
                   + '_'.join([options.labels[1], 'unique.bed']))
    (c-a-b).moveto(options.outputfolder + options.outputid
                   + '_'.join([options.labels[2], 'unique.bed']))
