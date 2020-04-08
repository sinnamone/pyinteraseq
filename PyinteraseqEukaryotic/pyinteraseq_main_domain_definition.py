import optparse
from pyinteraseq_domains_definition import *

parser = optparse.OptionParser(usage='python %prog pyinteraseq_main_domain_definition.py', version='1.0',)
input_opts = optparse.OptionGroup(
    parser, 'Input Options',
    'Input file.',
    )
input_opts.add_option('--mappingoutput', action="store", dest="mappingoutput", default=None,
                      help='BAM sorted generated from pyinteraseq_main_mapping.py')

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
reference_opts.add_option('--annotation', action="store", dest="annotation", default=None,
                          help='Annotation File(.gff|.bed)')
reference_opts.add_option('--genome', action="store", dest="genome", default=None,
                          help='Genome file')
reference_opts.add_option('--log', action="store", dest="log", default=None,
                          help='Path/filelog.log.')
parser.add_option_group(reference_opts)


reference_opts = optparse.OptionGroup(
    parser, 'Advanced Options',
    'Options for advanced analysis.',
    )
reference_opts.add_option('--threshold', action="store", dest="threshold", type="int",
                          default=1, help='Percentile threshold')
parser.add_option_group(reference_opts)

options, args = parser.parse_args()


if __name__ == '__main__':
    DictInfo = dict()
    DomainDefinitionClass = DomainsDefinition(optparseinstance=options)
    DictInfo["depthcoverage"] = DomainDefinitionClass.depthcoverage()
    DictInfo["breadthcoverage"] = DomainDefinitionClass.breadthcoverage()
    DictInfo["bam2bed"] = DomainDefinitionClass.bam2tabular()
    DictInfo["readcount"] = DomainDefinitionClass.bedtoolscoverage(DictInfo["bam2bed"])
    DictInfo["maxdepth"] = DomainDefinitionClass.groupbydepth(DictInfo["breadthcoverage"])
    DictInfo["percentile"] = DomainDefinitionClass.percentile(DictInfo["depthcoverage"])
    DictInfo["startenddef"] = DomainDefinitionClass.startenddefinition(DictInfo["percentile"])
    DictInfo["domainnotmerged"] = DomainDefinitionClass.domainparsing(DictInfo["startenddef"])
    DictInfo["domaindefinition"] = DomainDefinitionClass.parsingoutput(DictInfo["domainnotmerged"],
                                                                       DictInfo["readcount"])
    DomainDefinitionClass.cleantempfile()
