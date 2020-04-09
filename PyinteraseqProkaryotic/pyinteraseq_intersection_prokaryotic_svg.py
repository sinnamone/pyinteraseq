import optparse
from output_message import *
import pybedtools
import sys
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib_venn import venn3
import pandas as pd
import traceback
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
interection = 'intersection.txt'
unique = 'unique.txt'
header = ["#Chr", "CloneStart", "CloneEnd", "CloneLength", "Start", "End", "GeneID",
          "logFC", "PValue", "AdjPValue", "Strand", "Description", "NuclSeq"]
for i in range(len(options.selections)):
    count += 1
if options.log is None:
    filelog = options.outputfolder + options.outputid + '_intersection.log'
elif options.log is not None:
    filelog = options.log

if count == 0:
    filelogopen = open(str(filelog), 'a', 0)
    filelogopen.write(msg85)
    sys.exit(1)
elif count == 1:
    filelogopen = open(str(filelog), 'a', 0)
    filelogopen.write(msg86)
    sys.exit(1)
elif count == 2:
    try:
        filelogopen = open(str(filelog), 'a', 0)
        a = pybedtools.BedTool(options.selections[0])
        b = pybedtools.BedTool(options.selections[1])
        # intersection between two list
        aintersectb = a.intersect(b, u=True)
        # unique for list A
        aunique = a.intersect(b, v=True)
        # unique for list B
        bunique = b.intersect(a, v=True)
	# Giorgio
        # v = venn2(subsets=((a - b).count(), (a + b).count(), (b - a).count()),
                  # set_labels=(options.labels[0], options.labels[1]))
        v = venn2(subsets=((a - b).count(), (a + b).count(), (b - a).count()),
                  set_labels=('A', 'B'))
        plt.savefig(options.outputfolder + options.outputid + '.svg', format='svg')
        plt.savefig(options.outputfolder + options.outputid + '.png', format='png')
        # write files
        dfa_b = pd.read_table((a - b).fn, header=None)
        dfa_b[3] = dfa_b[12]
        dfa_b[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]].to_csv(options.outputfolder + options.outputid + '_' +
                                                                 '_'.join([options.labels[0], unique]),
                                                                 sep="\t", header=header, index=False)
        #
        dfb_a = pd.read_table((b - a).fn, header=None)
        dfb_a[3] = dfb_a[12]
        dfb_a[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]].to_csv(options.outputfolder + options.outputid + '_' +
                                                                 '_'.join([options.labels[0], options.labels[1], interection]),
                                                                 sep="\t", header=header, index=False)
        #
        dfAB = pd.read_table((a + b).fn, header=None)
        dfAB[3] = dfAB[12]
        dfAB[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]].to_csv(options.outputfolder + options.outputid + '_' +
                                                                '_'.join( [options.labels[1], unique]),
                                                                sep="\t", header=header, index=False)

        # svg file for venn2
        svgfile = options.outputfolder + options.outputid + '.svg'

        # read the content
        svg_read = open(svgfile, "r")
        svg_content = svg_read.read()
        svg_read.close()

        # ng-click actions
        actions = ['Component.Ctrl.GetOutputFile(' + str(index) +')' for index in range(3)]

        # insert the actions
        svg_content = svg_content.replace('<svg',   '<svg class="center-block img-responsive"')

        svg_content = svg_content.replace('<g id="patch_2">',   '<g id="patch_2" style="cursor: pointer" ng-click="'+ actions[0] +'">')
        svg_content = svg_content.replace('<g id="patch_3">',   '<g id="patch_3" style="cursor: pointer" ng-click="'+ actions[1] +'">')
        svg_content = svg_content.replace('<g id="patch_4">',   '<g id="patch_4" style="cursor: pointer" ng-click="'+ actions[2] +'">')

        # open and rewrite content
        svgfile = options.outputfolder + options.outputid + '.js_svg'

        svg_write = open(svgfile, "w")
        svg_write.write(svg_content)
        svg_read.close()
    except traceback:
        filelogopen = open(str(filelog), 'a', 0)
        filelogopen.write(msg155)
        sys.exit(1)
    else:
        filelogopen.write(msg156)

elif count == 3:
    try:
        filelogopen = open(str(filelog), 'a', 0)
        a = pybedtools.BedTool(options.selections[0])
        b = pybedtools.BedTool(options.selections[1])
        c = pybedtools.BedTool(options.selections[2])
	# GIORGIO
        # v = venn3(subsets=(
            # (a - b - c).count(), (b - a - c).count(), (a + b - c).count(), (c - a - b).count(), (c + a - b).count(),
            # (c + b - a).count(), (a + b + c).count()), set_labels=(options.labels[0], options.labels[1], options.labels[2]))
        v = venn3(subsets=(
            (a - b - c).count(), (b - a - c).count(), (a + b - c).count(), (c - a - b).count(), (c + a - b).count(),
            (c + b - a).count(), (a + b + c).count()), set_labels=('A', 'B', 'C'))
        plt.savefig(options.outputfolder + options.outputid + '.svg', format='svg')
        plt.savefig(options.outputfolder + options.outputid + '.png', format='png')
        #
        dfAB = pd.read_table((a + b).fn, header=None)

        dfAB.to_csv(options.outputfolder + options.outputid + '_'
                                                                + '_'.join([options.labels[0], options.labels[1], interection]), sep="\t",
                                                                header=header, index=False)
        #
        dfAC = pd.read_table((a + c).fn, header=None)

        dfAC.to_csv(options.outputfolder + options.outputid + '_' +
                                                                '_'.join([options.labels[0], options.labels[2],
                                                                          interection]), sep="\t",
                                                                header=header,  index=False)
        #
        dfBC = pd.read_table((b + c).fn, header=None)

        dfBC.to_csv(options.outputfolder + options.outputid + '_' +
                                                                '_'.join([options.labels[1], options.labels[2], interection]), sep="\t",
                                                                header=header, index=False)
        #
        dfABC = pd.read_table((a + b + c).fn, header=None)

        dfABC.to_csv(options.outputfolder + options.outputid + '_' +
                                                                 '_'.join([options.labels[0], options.labels[1], options.labels[2] + '_' + interection]), sep="\t",
                                                                 header=header, index=False)
        #
        dfa_b_c = pd.read_table((a - b - c).fn, header=None)

        dfa_b_c.to_csv(options.outputfolder + options.outputid + '_' +
                                                                   '_'.join([options.labels[0], unique]), sep="\t",
                                                                   header=header, index=False)
        #
        dfb_a_c = pd.read_table((b - a - c).fn, header=None)

        dfb_a_c.to_csv(options.outputfolder + options.outputid + '_' +
                                                                   '_'.join([options.labels[1], unique]), sep="\t",
                                                                   header=header, index=False)
        #
        dfc_a_b = pd.read_table((c - a - b).fn,  header=None)

        dfc_a_b.to_csv(options.outputfolder + options.outputid + '_' +
                                                                   '_'.join([options.labels[2], unique]), sep="\t",
                                                                   header=header, index=False)

        # svg file for venn3
        svgfile = options.outputfolder + options.outputid + '.svg'

        # read the content
        svg_read = open(svgfile, "r")
        svg_content = svg_read.read()
        svg_read.close()

        # ng-click actions
        actions = ['Component.Ctrl.GetOutputFile(' + str(index) +')' for index in range(7)]

        # insert the actions
        svg_content = svg_content.replace('<svg',   '<svg class="center-block img-responsive"')

        svg_content = svg_content.replace('<g id="patch_2">',   '<g id="patch_2" style="cursor: pointer" ng-click="'+ actions[0] +'">')
        svg_content = svg_content.replace('<g id="patch_3">',   '<g id="patch_3" style="cursor: pointer" ng-click="'+ actions[1] +'">')
        svg_content = svg_content.replace('<g id="patch_4">',   '<g id="patch_4" style="cursor: pointer" ng-click="'+ actions[3] +'">')
        svg_content = svg_content.replace('<g id="patch_5">',   '<g id="patch_5" style="cursor: pointer" ng-click="'+ actions[2] +'">')
        svg_content = svg_content.replace('<g id="patch_6">',   '<g id="patch_6" style="cursor: pointer" ng-click="'+ actions[4] +'">')
        svg_content = svg_content.replace('<g id="patch_7">',   '<g id="patch_7" style="cursor: pointer" ng-click="'+ actions[5] +'">')
        svg_content = svg_content.replace('<g id="patch_8">',   '<g id="patch_8" style="cursor: pointer" ng-click="'+ actions[6] +'">')

        # open and rewrite content
        svgfile = options.outputfolder + options.outputid + '.js_svg'

        svg_write = open(svgfile, "w")
        svg_write.write(svg_content)
        svg_read.close()
    except traceback:
        filelogopen = open(str(filelog), 'a', 0)
        filelogopen.write(msg155)
        sys.exit(1)
    else:
        filelogopen.write(msg156)

