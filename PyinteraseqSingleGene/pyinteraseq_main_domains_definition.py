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
parser.add_option_group(reference_opts)

options, args = parser.parse_args()

if __name__ == '__main__':
    DictInfo = dict()
    # 
    DomDef = DomainsDefinition(optparseinstance=options)
    DomDef.inputinformationappen()
    # count lenght gene
    DictInfo["fastalength"] = DomDef.contafasta(ref=options.fastasequence)
    # creation of genome file
    DictInfo["genomefile"] = DomDef.genomefile()
    # convert balstn to bed
    DictInfo["backgroundbed"] = DomDef.blastntobed(
        blastnoutput=options.backgroundmappingoutput, idex='background')
    DictInfo["targetbed"] = DomDef.blastntobed(
        blastnoutput=options.targetmappingoutput, idex='target')
    # conversion bed to bam
    DictInfo["backgroundbam"] = DomDef.bedtobam(
        bedinp=DictInfo["backgroundbed"], genomefile=DictInfo["genomefile"], idex="background")
    DictInfo["targetbam"] = DomDef.bedtobam(
        bedinp=DictInfo["targetbed"], genomefile=DictInfo["genomefile"], idex="target")
    #
    DictInfo["backgroundbamaveragedepth"] = DomDef.averagedepth(
        bam=DictInfo["backgroundbam"], lung=DictInfo["fastalength"])
    DictInfo["targetbamaveragedepth"] = DomDef.averagedepth(
        bam=DictInfo["targetbam"], lung=DictInfo["fastalength"])
    #
    DictInfo["downsampledbam"] = DomDef.downsampling(
        depb=DictInfo["backgroundbamaveragedepth"],
        dept=DictInfo["targetbamaveragedepth"],
        bamb=DictInfo["backgroundbam"],
        bamt=DictInfo["targetbam"])
    #
    if DictInfo["downsampledbam"][0] == "background":
        DictInfo["coveragebackground"] = DomDef.depthcount(
            bam=DictInfo["downsampledbam"][1], prefix=DictInfo["downsampledbam"][0])
        DictInfo["coveragetarget"] = DomDef.depthcount(
            bam=DictInfo["targetbam"], prefix="target")
    else:
        DictInfo["coveragetarget"] = DomDef.depthcount(
            bam=DictInfo["downsampledbam"][1], prefix=DictInfo["downsampledbam"][0])
        DictInfo["coveragebackground"] = DomDef.depthcount(
            bam=DictInfo["backgroundbam"], prefix="background")
    DictInfo["intervalsdata"] = DomDef.intervaldata(
        bedb=DictInfo["coveragebackground"], bedt=DictInfo["coveragetarget"])
    DictInfo["transpose"] = DomDef.trasposedomains(
        intervals=DictInfo["intervalsdata"])
    DictInfo["filtered"] = DomDef.filtering_domain(DictInfo["transpose"])
    for i in "targetbed", "targetbam", "coveragetarget", "backgroundbam", "coveragebackground", "intervalsdata", "backgroundbed", "genomefile", "transpose":
        os.remove(DictInfo[i])
    os.remove(DictInfo["downsampledbam"][1])


