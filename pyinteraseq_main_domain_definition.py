import optparse
from pyinteraseq_mapping import *
from pytinteraseq_domains_definition import *

parser = optparse.OptionParser(usage='python %prog pyinteraseq_main_definition.py', version='1.0',)
input_opts = optparse.OptionGroup(
    parser, 'Input Options',
    'Input file.',
    )
input_opts.add_option('--mappingoutput', action="store", dest="mappingoutput", default=None,
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
query_opts.add_option('--log', action="store", dest="log", default=None,
                      help='Log file.')
parser.add_option_group(query_opts)

reference_opts = optparse.OptionGroup(
    parser, 'Reference Options',
    'Options for reference FastA and Annotation.',
    )
reference_opts.add_option('--fastasequence', action="store", dest="fastasequence", default=None,
                          help='Genome sequence fasta file.(.fasta|.fna|.fa)')
reference_opts.add_option('--annotation', action="store", dest="annotation", default=None,
                          help='Annotation File(.gff|.bed)')
reference_opts.add_option('--chromosomename', action="store", dest="chromosomename", default=None,
                          help='Chromosome Name.(NC_XXXX)')
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
    DictFile = dict()
    outformat6 = '6 sseqid sstart send qseqid score sstrand'
    # Conversion Filter Blastn Table in Fasta
    DictInfo["blastoutputnohashfilteredfasta"] = DomainsDefinition(optparseinstance=options).tab2fasta(
        tabular=options.mappingoutput, prefixoutput="_blastnfiltered")
    # Clustering steps calling script Pick_otus
    DictInfo["clustering"] = DomainsDefinition(optparseinstance=options).clustering(
        blastnout=DictInfo["blastoutputnohashfilteredfasta"], prefixoutput="_blastnfiltered")
    # Pick most representative sequence for each cluster
    DictInfo["pickedreads"] = DomainsDefinition(optparseinstance=options).pickrepseq(
        pickotus=DictInfo["clustering"], fasta=DictInfo["blastoutputnohashfilteredfasta"])
    # Sed function
    DictInfo["pickedreadscleand"] = DomainsDefinition(optparseinstance=options).pysed(
        DictInfo["pickedreads"], '_clean.fasta', '-', '')
    # Mapping most representative clone against genome to identify which gene are intersted
    DictInfo["blastedclones"] = DomainsDefinition(optparseinstance=options).callmultiblastn(
        fasta=options.fastasequence,
        multifasta=DictInfo["pickedreadscleand"],
        outputformat=outformat6,
        suffix='_blastnclones.tab')
    # Parser that take as input blastn format 6 and create a standard bed6
    DictInfo["bedparsed"] = DomainsDefinition(optparseinstance=options).bedparsing(DictInfo["blastedclones"])
    # Filtering the output of Bedtools annotate (call inside the function) using a flot percentage of overlap
    DictInfo["clonesannotated"] = DomainsDefinition(optparseinstance=options).bedtoolsannotatefiltering(
        # Bedtools annotate to identify the clones inside the CDS
        DomainsDefinition(optparseinstance=options).bedtoolsannotate(
            DictInfo["bedparsed"], options.annotation), options.overlapintersect)
    # Create the count file for each intervals
    DictInfo["clustercount"] = DomainsDefinition(optparseinstance=options).clonescount(
        DictInfo["clustering"])
    # Merge BED6 with count table
    DictInfo["clonescounted"] = DomainsDefinition(optparseinstance=options).mergingcount(
        DictInfo["bedparsed"], DictInfo["clustercount"])
    #  Filtering the domains that are covered less than 10 clones
    DictInfo["clonescountedfiltered"] = DomainsDefinition(optparseinstance=options).filteringclonescount(
        DictInfo["clonescounted"], 10)
    # clone merge to get domains
    DictInfo["clonescountedmerged"] = DomainsDefinition(optparseinstance=options).pybedtoolsmerge(
        DictInfo["clonescountedfiltered"])
    # Get fasta from intervels
    DictInfo["clonesmergedfasta"] = DomainsDefinition(optparseinstance=options).pybedtoolstofasta(
        pybedtoolsmergeoutput=DictInfo["clonescountedmerged"], fastqsequence=options.fastasequence)
    # Add description present in ptt 0.7 is the overlap of intersect cds
    DictInfo["tabwithdescription"] = DomainsDefinition(optparseinstance=options).adddescription(
        clonesmerged=DictInfo["clonescountedmerged"], annotation=options.annotation,
        percthr=options.overlapintersect)
    # Conversion Fasta ==>Tabular
    DictInfo["clonenseqfasta"] = DomainsDefinition(optparseinstance=options).fasta2tabular(
        imp=DictInfo["clonesmergedfasta"], prefix='_clonestabular')
    # Add sequence to output table
    DictFile["tabwithsequence"] = DomainsDefinition(optparseinstance=options).addsequence(
        outputfromdescription=DictInfo["tabwithdescription"], outputfasta2tab=DictInfo["clonenseqfasta"])
    # print DictInfo
    # DomainsDefinition(optparseinstance=options).cleantemporaryfilesinglend(DictInfo)
