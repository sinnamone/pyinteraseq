import optparse
from pyinteraseq_domains_definition import *

parser = optparse.OptionParser(usage='python %prog pyinteraseq_main_definition.py', version='1.0',)
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
                          help='Genome File(.gff|.bed)')
reference_opts.add_option('--log', action="store", dest="log", default=None,
                          help='Path/filelog.log.')
parser.add_option_group(reference_opts)


reference_opts = optparse.OptionGroup(
    parser, 'Advanced Options',
    'Options for advanced analysis.',
    )
reference_opts.add_option('--minclonelength', action="store", dest="minclonelength", default='100',
                          help='Minumum clones length.')
reference_opts.add_option('--overlapintersect', action="store", dest="overlapintersect", type="float",
                          default=0.7, help='Parameters -f of bedtools intersect.')
reference_opts.add_option('--opengap', action="store", dest="opengap", type="int",
                          default=1, help='Open-gap allowed.')
reference_opts.add_option('--mismatch', action="store", dest="mismatch", type="float",
                          default=3.0, help='Percentage of Mismatch allowed.')
reference_opts.add_option('--thread', action="store", dest="thread", default='1',
                          help='Number of thread.')
parser.add_option_group(reference_opts)

options, args = parser.parse_args()

if __name__ == '__main__':
    DictInfo = dict()
    outformat6 = '6 sseqid sstart send qseqid score sstrand'
    DomainDefinitionClass = DomainsDefinition(optparseinstance=options)
    DictInfo["bam2tabular"] = DomainDefinitionClass.bam2tabular()
    DictInfo["fasta2tabular"] = DomainDefinitionClass.seqrename(tabular=DictInfo["bam2tabular"])
    DictInfo["tabular2fasta"] = DomainDefinitionClass.tab2fasta(tabular=DictInfo["fasta2tabular"],
                                                                prefixoutput="_renamed")
    DictInfo["clustering"] = DomainDefinitionClass.clustering(bowtie2out=DictInfo["tabular2fasta"])
    # Pick most representative sequence for each cluster
    DictInfo["pickedreads"] = DomainDefinitionClass.pickrepseq(
        pickotus=DictInfo["clustering"], fasta=DictInfo["tabular2fasta"])
    # Sed function
    DictInfo["pickedreadscleand"] = DomainDefinitionClass.pysed(
        DictInfo["pickedreads"], '_clean.fasta', '-', '')
    # Mapping most representative clone against genome to identify which gene are intersted
    DictInfo["blastedclones"] = DomainDefinitionClass.callmultiblastn(
        fasta=DomainDefinitionClass.fastasequence,
        multifasta=DictInfo["pickedreadscleand"],
        outputformat=outformat6,
        suffix='_blastnclones.tab')
    # Parser that take as input blastn format 6 and create a standard bed6
    DictInfo["bedparsed"] = DomainDefinitionClass.bedparsing(DictInfo["blastedclones"])
    DictInfo["bedfixed"] = DomainDefinitionClass.bedannotateparsing(bedparsed=DictInfo["bedparsed"])
    DictInfo["cloneschecked"] = DomainDefinitionClass.bedtoolsannotate(
        clonesformatbed=DictInfo["bedfixed"],
        annotation=DomainDefinitionClass.annotation)
    # Bedtools annotate to identify the clones inside the CDS
    DictInfo["clonesannotated"] = DomainDefinitionClass.bedtoolsannotatefiltering(
        bedtoolsannotateout=DictInfo["cloneschecked"],
        percthr=DomainDefinitionClass.overlapintersect)
    #
    DictInfo["clustercount"] = DomainDefinitionClass.clonescount(
        DictInfo["clustering"])
    # Merge BED6 with count table
    DictInfo["clonescounted"] = DomainDefinitionClass.mergingcount(
        DictInfo["bedfixed"], DictInfo["clustercount"])
    # #  Filtering the domains that are covered less than 10 clones
    DictInfo["clonescountedfiltered"] = DomainDefinitionClass.filteringclonescount(
        DictInfo["clonescounted"], 5)
    # clone merge to get domains
    DictInfo["clonescountedmerged"] = DomainDefinitionClass.pybedtoolsmerge(
        DictInfo["clonescountedfiltered"])
    # Get fasta from intervels
    DictInfo["clonesmergedfasta"] = DomainDefinitionClass.pybedtoolstofasta(
        pybedtoolsmergeoutput=DictInfo["clonescountedmerged"], fastqsequence=options.fastasequence)
    # Add description present in ptt 0.7 is the overlap of intersect cds
    DictInfo["tabwithdescription"] = DomainDefinitionClass.adddescription(
        clonesmerged=DictInfo["clonescountedmerged"], annotation=options.annotation,
        percthr=options.overlapintersect)
    # Add sequence to output table
    DictInfo["tabwithsequence"] = DomainDefinitionClass.addsequence(
        outputfromdescription=DictInfo["tabwithdescription"], outputfasta2tab=DictInfo["clonesmergedfasta"])
    # DomainDefinitionClass.cleantempfile()
