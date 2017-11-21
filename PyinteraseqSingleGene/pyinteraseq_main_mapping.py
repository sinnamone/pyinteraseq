import optparse
from pyinteraseq_mapping import *

parser = optparse.OptionParser(usage='python %prog pyinteraseq_main_mapping.py', version='1.0',)
parser.add_option('--readforwardtrimmed', action="store", dest="readforwardtrimmed", default=None,
                  help='Read dataset input forward')
parser.add_option('--readreversetrimmed', action="store", dest="readreversetrimmed", default=None,
                  help='Read dataset input reverse',)
parser.add_option('--readforwardtrimmedtype', type='choice', choices=['fastq', 'fasta'], default=None,
                  help='Read dataset input forward')
parser.add_option('--readreversetrimmedtype', type='choice', choices=['fastq', 'fasta'], default=None,
                  help='Read dataset input reverse')
parser.add_option('--sampletype', type='choice', choices=['background', 'target'], default=None,
                  help='Select type of dataset')

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
    parser, 'Primers Options',
    'Options for primer sequence and trimming.',
    )
reference_opts.add_option('--fastasequence', action="store", dest="fastasequence", default=None,
                          help='Genome sequence fasta file.(.fasta|.fna|.fa)')
reference_opts.add_option('--genename', action="store", dest="genename", default=None,
                          help='String with gene name')
parser.add_option_group(reference_opts)

advance_opts = optparse.OptionGroup(
    parser, 'Advanced Options',
    'Options for advanced analysis.',
    )
advance_opts.add_option('--minclonelength', action="store", dest="minclonelength", default='100',
                        help='Minumum clones length.')
advance_opts.add_option('--overlapintersect', action="store", dest="overlapintersect", type="float",
                        default=0.7, help='Parameters -f of bedtools intersect.')
advance_opts.add_option('--opengap', action="store", dest="opengap", type="int",
                        default=1, help='Open-gap allowed.')
advance_opts.add_option('--mismatch', action="store", dest="mismatch", type="float",
                        default=3.0, help='Percentage of Mismatch allowed.')
parser.add_option_group(advance_opts)

options, args = parser.parse_args()

if __name__ == '__main__':
    DictInfo = dict()
    DictFile = dict()
    outformat7 = '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq'
    MergeFileList = []
    if options.readforwardtrimmedtype == 'fastq':
        # Conversion Fastq==>Fasta Paired-end
        if options.readreversetrimmed is not None:
            DictInfo["FastaReadsForward"] = BlastNlucleotide(optparseinstance=options).fastq2fasta(
                fastq=options.readforwardtrimmed, nameid="_forward")
            # Conversion Fastq<==>Fasta Paired-end
            DictInfo["FastaReadsReverse"] = BlastNlucleotide(optparseinstance=options).fastq2fasta(
                fastq=options.readreversetrimmed, nameid="_reverse")
        else:
            # Conversion Fastq==>Fasta Single-end
            DictInfo["FastaReadsForward"] = BlastNlucleotide(optparseinstance=options).fastq2fasta(
                fastq=options.readforwardtrimmed, nameid="_forward")
        # Conversion Fasta==>Tabular forward
        DictInfo["TabularReadsForward"] = BlastNlucleotide(optparseinstance=options).fasta2tabular(
            imp=DictInfo["FastaReadsForward"], prefix="forward")
        if options.readreversetrimmed is not None:
            # Conversion Fasta==>Tabular reverse
            DictInfo["TabularReadsReverse"] = BlastNlucleotide(optparseinstance=options).fasta2tabular(
                imp=DictInfo["FastaReadsReverse"], prefix="reverse")
    elif options.readforwardtrimmedtype == 'fasta':
        # Conversion FastA<==>Tabular Single-End
        DictInfo["TabularReadsForward"] = BlastNlucleotide(optparseinstance=options).fasta2tabular(
            imp=options.readforwardtrimmed, prefix="forward")
        if options.readreversetrimmed is not None:
            DictInfo["TabularReadsReverse"] = BlastNlucleotide(optparseinstance=options).fasta2tabular(
                imp=options.readreversetrimmed, prefix="reverse")
    # Rename Sequence ID forward
    DictInfo["TabularRenamedForward"] = BlastNlucleotide(optparseinstance=options).seqrename(
        tabular=DictInfo["TabularReadsForward"], readirection="forward")
    # Renanme Sequence ID reverse
    if options.readreversetrimmed is not None:
        DictInfo["TabularRenamedReverse"] = BlastNlucleotide(optparseinstance=options).seqrename(
            tabular=DictInfo["TabularReadsReverse"], readirection="reverse")
    DictInfo["FastaRenamedForward"] = BlastNlucleotide(optparseinstance=options).tab2fasta(
        DictInfo["TabularRenamedForward"], "_forward")
    if options.readreversetrimmed is not None:
        # Conversion Tabular <==>Fasta reverse
        DictInfo["FastaRenamedReverse"] = BlastNlucleotide(optparseinstance=options).tab2fasta(
            DictInfo["TabularRenamedReverse"], "_reverse")
        # Merging Pairs
        DictInfo["Trimmedreadconcatenated"] = BlastNlucleotide(optparseinstance=options).concatenateforrev(
            [DictInfo["FastaRenamedForward"], DictInfo["FastaRenamedReverse"]])
        # Mapping steps for Paired-End
        DictInfo["blastoutput"] = BlastNlucleotide(optparseinstance=options).callmultiblastn(
            fasta=options.fastasequence,
            multifasta=DictInfo["Trimmedreadconcatenated"],
            outputformat=outformat7,
            suffix='_blastn.txt')
    else:
        # Mapping steps for Single-End
        DictInfo["blastoutput"] = BlastNlucleotide(optparseinstance=options).callmultiblastn(
            fasta=options.fastasequence,
            multifasta=DictInfo["FastaRenamedForward"],
            outputformat=outformat7,
            suffix='_blastn.txt')
    # Clean output of hash
    DictInfo["blastoutputnohash"] = BlastNlucleotide(optparseinstance=options).hashclean(
        blastnout=DictInfo["blastoutput"], prefix="_blastn_nohash")
    # Filter reads steps (NO open-gaps, mismatch)
    DictInfo["blastoutputnohashfiltered"] = BlastNlucleotide(optparseinstance=options).blastnfiltering(
        blastnout=DictInfo["blastoutputnohash"])
