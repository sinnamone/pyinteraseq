import optparse
from output_message import *
import pybedtools
import sys
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

parser = optparse.OptionParser()
parser.add_option('-i', action="append", dest='inputs', default=[])
parser.add_option('-l', action="append", dest='labels', default=[])
parser.add_option('--outputfolder', action="store", dest="outputfolder", default=None,
                      help='Output folder.')
parser.add_option('--outputid', action="store", dest="outputid", default=None,
                      help='Output ID.')
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
    aintersectb = a.intersect(b, f=0.50, r=True, u=True)
    # unique for list A
    aunique = a.intersect(b, f=0.50, r=True, v=True)
    # unique for list B
    bunique = b.intersect(a, f=0.50, r=True, v=True)
    v = venn2(subsets=(len(aunique), len(aintersectb), len(bunique)), set_labels=(options.labels[0], options.labels[1]))
    v.get_patch_by_id('10').set_alpha(0.4)
    v.get_patch_by_id('01').set_alpha(1.0)
    v.get_patch_by_id('11').set_alpha(0.7)
    plt.savefig(options.outputfolder + options.outputid + '.pdf', format='pdf')
elif count == 3:
    a = pybedtools.BedTool(options.inputs[0])
    b = pybedtools.BedTool(options.inputs[1])
    c = pybedtools.BedTool(options.inputs[2])
    # intersection between A and B
    aintersectb = a.intersect(b, f=0.50, r=True, u=True)
    # intersection between A and C
    aintersectc = a.intersect(c, f=0.50, r=True, u=True)
    # intersection between B and C
    bintersectc = b.intersect(c, f=0.50, r=True, u=True)
    # intersection between AB e AC
    abintersectac = aintersectb.intersect(aintersectc, f=0.50, r=True, u=True)
    # common all
    commonall = abintersectac.intersect(bintersectc, f=0.50, r=True, u=True)
    # unique between A and B
    abunique = a.intersect(b, f=0.50, r=True, v=True)
    # unique list A
    aunique = abunique.intersect(c, f=0.50, r=True, v=True)
    # unique between B and C
    baunique = b.intersect(c, f=0.50, r=True, v=True)
    # unique list B
    bunique = baunique.intersect(a, f=0.50, r=True, v=True)
    # unique between C and B
    cbunique = c.intersect(b, f=0.50, r=True, v=True)
    # unique list C
    cunique = cbunique.intersect(a, f=0.50, r=True, v=True)
    print (len(aunique), len(abunique), len(bunique), len(cunique), len(aintersectc), len(bintersectc), len(commonall))
    v = venn2(subsets=(len(aunique), len(abunique), len(bunique), len(cunique), len(aintersectc),
                       len(bintersectc), len(commonall)), set_labels=(options.labels[0], options.labels[1], options.labels[2]))
    plt.savefig(options.outputfolder + options.outputid + '.pdf', format='pdf')
