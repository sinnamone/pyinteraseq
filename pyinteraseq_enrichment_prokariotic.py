from pyinteraseq_genomefileparsing import *
import optparse

parser = optparse.OptionParser(usage='python %prog Main PyInteract', version='1.0',)
parser.add_option('--readforwardgenomic', action="store", dest="readforwardgenomic", default=None,
                  help='Read dataset input forward')
parser.add_option('--readreversegenomic', action="store", dest="readreversegenomic", default=None,
                  help='Read dataset input reverse',)
parser.add_option('--readforwardtypegenomic', type='choice', choices=['fastq', 'fasta'], default=None,
                  help='Read dataset input forward')
parser.add_option('--readreversetypegenomic', type='choice', choices=['fastq', 'fasta'], default=None,
                  help='Read dataset input reverse')
parser.add_option('--readforwardtarget', action="store", dest="readforwardtarget", default=None,
                  help='Read dataset input forward')
parser.add_option('--readreversegenomictarget', action="store", dest="readreversetarget", default=None,
                  help='Read dataset input reverse',)
parser.add_option('--readforwardtypegenomictarget', type='choice', choices=['fastq', 'fasta'], default=None,
                  help='Read dataset input forward')
parser.add_option('--readreversetypegenomictarget', type='choice', choices=['fastq', 'fasta'], default=None,
                  help='Read dataset input reverse')
parser.add_option('--domaindetectiontarget', action="store", dest="readreversetarget", default=None,
                  help='Output from pyinterseq_domains_definition')

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
    parser, 'Primers Options',
    'Options for primer sequence and trimming.',
    )
reference_opts.add_option('--fastasequence', action="store", dest="fastasequence", default=None,
                          help='Genome sequence fasta file.')
reference_opts.add_option('--annotation', action="store", dest="annotation", default=None,
                          help='Annotation File')
parser.add_option_group(reference_opts)

parser.add_option_group(query_opts)

reference_opts = optparse.OptionGroup(
    parser, 'Advanced Options',
    'Options for advanced analysis.',
    )
reference_opts.add_option('--thread', action="store", dest="thread", default='2',
                          help='Number of thread.')
parser.add_option_group(reference_opts)

options, args = parser.parse_args()

if __name__ == '__main__':
    DictInfo = dict()
    DictInfo["fasta"] = GenomeFile(optparseinstance=options).fastareference()
    #
    DictInfo["annotation"] = AnnotationFile(optparseinstance=options).annotationbuild()
    #
