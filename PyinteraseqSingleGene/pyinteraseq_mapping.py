import optparse
from pyinteraseq_mapping import *
from pyinteraseq_trimming import *
from pyinteraseq_genomefileparsing import *


parser = optparse.OptionParser(usage='python %prog pyinteraseq_main_mapping.py', version='1.0',)
parser.add_option('--readforward', action="store", dest="readforward", default=None,
                  help='Read dataset input forward')
parser.add_option('--readreverse', action="store", dest="readreverse", default=None,
                  help='Read dataset input reverse')

primer_opts = optparse.OptionGroup(
    parser, 'Primers Options',
    'Options for primer sequence and trimming.',
    )
primer_opts.add_option('--primer5forward', action="store", dest="primer5forward", default=None,
                       help='Output folder.')
primer_opts.add_option('--primer3forward', action="store", dest="primer3forward", default=None,
                       help='Output ID.')
primer_opts.add_option('--primer5reverse', action="store", dest="primer5reverse", default=None,
                       help='Output folder.')
primer_opts.add_option('--primer3reverse', action="store", dest="primer3reverse", default=None,
                       help='Output ID.')
parser.add_option_group(primer_opts)

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
                          help='Genome sequence fasta file.(.fasta|.fna|.fa)')
reference_opts.add_option('--thread', action="store", dest="thread", default='1',
                          help='Number of thread.')
reference_opts.add_option('--log', action="store", dest="log", default=None,
                          help='Number of thread.')
parser.add_option_group(reference_opts)

advance_opts = optparse.OptionGroup(
    parser, 'Advanced Options',
    'Options for advanced analysis.',
    )
advance_opts.add_option('--minclonelength', action="store", dest="minclonelength", default='100',
                        help='Minumum clones length.')
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
    InpClass = InputCheck(optparseinstance=options)
    DictInfo.update({
        "CutadaptPath": InpClass.cutadaptcheck(),
        "LogInfoAppended": InpClass.inputinformationappen()
    })
    # Parsing input fasta file
    DictInfo["fasta"] = GenomeFile(optparseinstance=options).fastareference()
    # Trimming steps for fastq input
    TrimmingSingleClass = TrimmingSingle(optparseinstance=options)
    TrimmingPairedClass = TrimmingPaired(optparseinstance=options)
    MappingClass = BlastNlucleotide(optparseinstance=options)
    readtype = InpClass.fastatesting(InpClass.readforward)
    seqtype = InpClass.sequencingtype
    if readtype == 'fastq':
        # Trimming paired end fastq
        if seqtype == "Paired-End":
            DictInfo["Trimmed5paired"] = TrimmingPairedClass.trimming5paired()
            DictInfo["FastaReadsForward"] = MappingClass.fastq2fasta(
                fastq=DictInfo["Trimmed5paired"][0], nameid="_forward")
            DictInfo["FastaReadsReverse"] = MappingClass.fastq2fasta(
                fastq=DictInfo["Trimmed5paired"][1], nameid="_reverse")
        elif seqtype == "Single-End":
            # Trimming Single-End fastq
            DictInfo["Trimmed5single"] = TrimmingSingleClass.trimming5single()
            DictInfo["FastaReadsForward"] = MappingClass.fastq2fasta(
                fastq=InpClass.readforward, nameid="_forward")
        else:
            sys.exit(1)
        # Conversion Fasta==>Tabular forward
        DictInfo["TabularReadsForward"] = MappingClass.fasta2tabular(
            imp=DictInfo["FastaReadsForward"], prefix="_forward")
        if seqtype == "Paired-End":
            # Conversion Fasta==>Tabular reverse
            DictInfo["TabularReadsReverse"] = MappingClass.fasta2tabular(
                imp=DictInfo["FastaReadsReverse"], prefix="_reverse")
    elif readtype == 'fasta':
        # Trimming paired end fasta
        if seqtype == "Paired-End":
            DictInfo["Trimmed5paired"] = TrimmingPairedClass.trimming5paired()
            DictInfo["TabularReadsForward"] = MappingClass.fasta2tabular(
                imp=InpClass.readforward, prefix="_forward")
            if options.readreverse is not None:
                DictInfo["TabularReadsReverse"] = MappingClass.fasta2tabular(
                    imp=InpClass.readreverse, prefix="_reverse")
        elif seqtype == "Single-End":
            # Trimming Single-End fasta
            DictInfo["Trimmed5single"] = TrimmingSingleClass.trimming5single()
            DictInfo["TabularReadsForward"] = MappingClass.fasta2tabular(
                imp=InpClass.readforward, prefix="_forward")
        else:
            sys.exit(1)
    # Rename Sequence ID forward
    DictInfo["TabularRenamedForward"] = MappingClass.seqrename(
        tabular=DictInfo["TabularReadsForward"], readirection="forward")
    # Renanme Sequence ID reverse
    if seqtype == "Paired-End":
        DictInfo["TabularRenamedReverse"] = MappingClass.seqrename(
            tabular=DictInfo["TabularReadsReverse"], readirection="reverse")
    DictInfo["FastaRenamedForward"] = MappingClass.tab2fasta(
        DictInfo["TabularRenamedForward"], "_forward")
    if seqtype == "Paired-End":
        # Conversion Tabular <==>Fasta reverse
        DictInfo["FastaRenamedReverse"] = MappingClass.tab2fasta(
            DictInfo["TabularRenamedReverse"], "_reverse")
        # Merging Pairs
        DictInfo["Trimmedreadconcatenated"] = MappingClass.concatenateforrev(
            [DictInfo["FastaRenamedForward"], DictInfo["FastaRenamedReverse"]])
        # Mapping steps for Paired-End
        DictInfo["blastoutput"] = MappingClass.callmultiblastn(
            fasta=options.fastasequence,
            multifasta=DictInfo["Trimmedreadconcatenated"],
            outputformat=outformat7,
            suffix='_blastn.txt')
    elif seqtype == "Single-End":
        # Mapping steps for Single-End
        DictInfo["blastoutput"] = MappingClass.callmultiblastn(
            fasta=options.fastasequence,
            multifasta=DictInfo["FastaRenamedForward"],
            outputformat=outformat7,
            suffix='_blastn.txt')
    else:
        sys.exit(1)
    # Clean output of hash
    DictInfo["blastoutputnohash"] = MappingClass.hashclean(
        blastnout=DictInfo["blastoutput"], prefix="_blastn_nohash")
    # Filter reads steps (NO open-gaps, mismatch)
    DictInfo["blastoutputnohashfiltered"] = MappingClass.blastnfiltering(
        blastnout=DictInfo["blastoutputnohash"])
    MappingClass.cleantempfile()
