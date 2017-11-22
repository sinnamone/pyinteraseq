import optparse
from pyinteraseq_mapping import *
from pyinteraseq_domains_definition import *

parser = optparse.OptionParser(usage='python %prog pyinteraseq_main_definition.py', version='1.0',)
input_opts = optparse.OptionGroup(
    parser, 'Input Options',
    'Input file.',
    )
input_opts.add_option('--backgroundmappingoutput', action="store", dest="backgroundmappingoutput", default=None,
                      help='Blastn output')
input_opts.add_option('--targetmappingoutput', action="store", dest="targetmappingoutput", default=None,
                      help='Blastn output')

parser.add_option_group(input_opts)

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
    parser, 'Reference Options',
    'Options for reference FastA and Annotation.',
    )
reference_opts.add_option('--fastasequence', action="store", dest="fastasequence", default=None,
                          help='Genome sequence fasta file.(.fasta|.fna|.fa)')
reference_opts.add_option('--genename', action="store", dest="genename", default=None,
                          help='String with gene name')
parser.add_option_group(reference_opts)

options, args = parser.parse_args()

if __name__ == '__main__':
    DictInfo = dict()
    DomainsDefinition(optparseinstance=options).inputinformationappen()
    # count lenght gene
    DictInfo["fastalength"] = DomainsDefinition(optparseinstance=options).contafasta(ref=options.fastasequence)
    # creation of genome file
    DictInfo["genomefile"] = DomainsDefinition(optparseinstance=options).genomefile()
    # convert balstn to bed
    DictInfo["backgroundbed"] = DomainsDefinition(optparseinstance=options).blastntobed(
        blastnoutput=options.backgroundmappingoutput, idex='background')
    DictInfo["targetbed"] = DomainsDefinition(optparseinstance=options).blastntobed(
        blastnoutput=options.targetmappingoutput, idex='target')
    # conversion bed to bam
    DictInfo["backgroundbam"] = DomainsDefinition(optparseinstance=options).bedtobam(
        bedinp=DictInfo["backgroundbed"], genomefile=DictInfo["genomefile"], idex="background")
    DictInfo["targetbam"] = DomainsDefinition(optparseinstance=options).bedtobam(
        bedinp=DictInfo["targetbed"], genomefile=DictInfo["genomefile"], idex="target")
    #
    DictInfo["backgroundbamaveragedepth"] = DomainsDefinition(optparseinstance=options).averagedepth(
        bam=DictInfo["backgroundbam"], lung=DictInfo["fastalength"])
    DictInfo["targetbamaveragedepth"] = DomainsDefinition(optparseinstance=options).averagedepth(
        bam=DictInfo["targetbam"], lung=DictInfo["fastalength"])
    #
    DictInfo["downsampledbam"] = DomainsDefinition(optparseinstance=options).downsampling(
        depb=DictInfo["backgroundbamaveragedepth"],
        dept=DictInfo["targetbamaveragedepth"],
        bamb=DictInfo["backgroundbam"],
        bamt=DictInfo["targetbam"])
    #
    if DictInfo["downsampledbam"][0] == "background":
        DictInfo["coveragebackground"] = DomainsDefinition(optparseinstance=options).depthcount(
            bam=DictInfo["downsampledbam"][1], prefix=DictInfo["downsampledbam"][0])
        DictInfo["coveragetarget"] = DomainsDefinition(optparseinstance=options).depthcount(
            bam=DictInfo["targetbam"], prefix="target")
    else:
        DictInfo["coveragetarget"] = DomainsDefinition(optparseinstance=options).depthcount(
            bam=DictInfo["downsampledbam"][1], prefix=DictInfo["downsampledbam"][0])
        DictInfo["coveragebackground"] = DomainsDefinition(optparseinstance=options).depthcount(
            bam=DictInfo["backgroundbam"], prefix="background")
    DomainsDefinition(optparseinstance=options).intervaldata(
        bedb=DictInfo["coveragebackground"], bedt=DictInfo["coveragetarget"])


