import optparse
from pyinteraseq_mapping import *
from pyinteraseq_trimming import *
from pyinteraseq_genomefileparsing import *
from output_message import *

parser = optparse.OptionParser(
    usage='python %prog pyinteraseq_main_mapping.py', version='1.0',)
parser.add_option('--readforward', action="store", dest="readforward", default=None,
                  help='Read dataset input forward')
parser.add_option('--readreverse', action="store", dest="readreverse", default=None,
                  help='Read dataset input reverse',)

query_opts = optparse.OptionGroup(
    parser, 'Output Options',
    'Options for the output destionation and name.',
)
query_opts.add_option('--outputfolder', action="store", dest="outputfolder", default=None,
                      help='Output folder.')
query_opts.add_option('--outputid', action="store", dest="outputid", default=None,
                      help='Output ID.')
parser.add_option_group(query_opts)

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

reference_opts = optparse.OptionGroup(
    parser, 'Primers Options',
    'Options for primer sequence and trimming.',
)
reference_opts.add_option('--fastasequence', action="store", dest="fastasequence", default=None,
                          help='Genome sequence fasta file.(.fasta|.fna|.fa)')
reference_opts.add_option('--thread', action="store", dest="thread", default='2',
                          help='Number of thread.')
reference_opts.add_option('--log', action="store", dest="log", default=None,
                          help='Log file.')
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
    InpClass = InputCheck(optparseinstance=options)
    MappingClass = BlastNlucleotide(optparseinstance=options)
    TrimSingle = TrimmingSingle(optparseinstance=options)
    TrimPaired = TrimmingPaired(optparseinstance=options)
    # input check
    if InpClass.readreverse is not None:
        if (InpClass.readforwardtype is "fastq") and (InpClass.readreversetype is "fastq"):
            DictInfo.update({
                "CutadaptPath": InpClass.cutadaptcheck(),
                "LogInfoAppended": InpClass.inputinformationappen()
            })
        elif (InpClass.readforwardtype is "fasta") and (InpClass.readreversetype is "fasta"):
            DictInfo.update({
                "CutadaptPath": InpClass.cutadaptcheck(),
                "LogInfoAppended": InpClass.inputinformationappen()
            })
        elif InpClass.readforwardtype is "fastq" and InpClass.readreversetype is "fasta":
            log = open(InpClass.inputfilelog, "a", 0)
            log.write(msg5)
        elif InpClass.readforwardtype is "fasta" and InpClass.readreversetype is "fastq":
            log = open(InpClass.inputfilelog, "a", 0)
            log.write(msg6)
    elif InpClass.readreverse is None:
        DictInfo.update({
            "CutadaptPath": InpClass.cutadaptcheck(),
            "LogInfoAppended": InpClass.inputinformationappen()
        })
    # start analysis
    log = open(InpClass.inputfilelog, "a", 0)
    if (InpClass.sequencingtype == "Single-End") and (InpClass.readforwardtype == "fastq"):
        if (InpClass.primer5forward is not None) and (InpClass.primer3forward is not None):
            # Trimming Single-End fastQ
            DictInfo["Trimmed5single"] = TrimSingle.trimming5single()
            DictInfo["FastaReadsForward"] = MappingClass.fastq2fasta(
                fastq=DictInfo["Trimmed5single"], nameid="_forward")
        elif (InpClass.primer5forward is None) and (InpClass.primer3forward is None):
            DictInfo["FastaReadsForward"] = MappingClass.fastq2fasta(
                fastq=InpClass.readforward, nameid="_forward")
            # Conversion Fasta==>Tabular forward
        else:
            log.write(msg116)
            sys.exit(1)
        DictInfo["TabularReadsForward"] = MappingClass.fasta2tabular(
            imp=DictInfo["FastaReadsForward"], prefix="forward")
        # Rename Sequence ID forward
        DictInfo["TabularRenamedForward"] = MappingClass.seqrename(
            tabular=DictInfo["TabularReadsForward"], readirection="forward")
        DictInfo["FastaRenamedForward"] = MappingClass.tab2fasta(
            DictInfo["TabularRenamedForward"], "_forward")
        # Mapping steps for Single-End
        DictInfo["blastoutput"] = MappingClass.callmultiblastn(
            fasta=options.fastasequence,
            multifasta=DictInfo["FastaRenamedForward"],
            outputformat=outformat7,
            suffix='_blastn.txt')
        # End Single fastq
    elif (InpClass.sequencingtype == "Single-End") and (InpClass.readforwardtype == "fasta"):
        if (InpClass.primer5forward is not None) and (InpClass.primer3forward is not None):
            # Trimming Single-End fastA
            DictInfo["Trimmed5single"] = TrimSingle.trimming5single()
        elif (InpClass.primer5forward is None) and (InpClass.primer3forward is None):
            DictInfo["Trimmed5single"] = InpClass.readforward
        else:
            log.write(msg116)
            sys.exit(1)
        DictInfo["TabularReadsForward"] = MappingClass.fasta2tabular(
            imp=DictInfo["Trimmed5single"], prefix="forward")
        # Rename Sequence ID forward
        DictInfo["TabularRenamedForward"] = MappingClass.seqrename(
            tabular=DictInfo["TabularReadsForward"], readirection="forward")
        DictInfo["FastaRenamedForward"] = MappingClass.tab2fasta(
            DictInfo["TabularRenamedForward"], "_forward")
        # Mapping steps for Single-End
        DictInfo["blastoutput"] = MappingClass.callmultiblastn(
            fasta=options.fastasequence,
            multifasta=DictInfo["FastaRenamedForward"],
            outputformat=outformat7,
            suffix='_blastn.txt')
        # End Single fastA
    elif (InpClass.sequencingtype == "Paired-End") and (InpClass.readforwardtype == "fastq"):
        if (InpClass.primer5forward is not None) and (InpClass.primer3forward is not None) and (InpClass.primer5reverse is not None) and (InpClass.primer3reverse is not None):
            DictInfo["Trimmed5paired"] = TrimPaired.trimming5paired()
            # Conversion Fastq<==>Fasta Paired-end
            DictInfo["FastaReadsForward"] = MappingClass.fastq2fasta(
                fastq=DictInfo["Trimmed5paired"][0], nameid="_forward")
            DictInfo["FastaReadsReverse"] = MappingClass.fastq2fasta(
                fastq=DictInfo["Trimmed5paired"][1], nameid="_reverse")
        elif (InpClass.primer5forward is None) and (InpClass.primer3forward is None) and (InpClass.primer5reverse is None) and (InpClass.primer3reverse is None):
            # Conversion Fastq<==>Fasta Paired-end
            DictInfo["FastaReadsForward"] = MappingClass.fastq2fasta(
                fastq=InpClass.readforward, nameid="_forward")
            DictInfo["FastaReadsReverse"] = MappingClass.fastq2fasta(
                fastq=InpClass.readreverse, nameid="_reverse")
        else:
            log.write(msg116)
            sys.exit(1)
        # Conversion  Fasta<==>Tabular
        DictInfo["TabularReadsForward"] = MappingClass.fasta2tabular(
            imp=DictInfo["FastaReadsForward"], prefix="forward")
        DictInfo["TabularReadsReverse"] = MappingClass.fasta2tabular(
            imp=DictInfo["FastaReadsReverse"], prefix="reverse")
        # Rename Sequence ID forward
        DictInfo["TabularRenamedForward"] = MappingClass.seqrename(
            tabular=DictInfo["TabularReadsForward"], readirection="forward")
        DictInfo["TabularRenamedReverse"] = MappingClass.seqrename(
            tabular=DictInfo["TabularReadsReverse"], readirection="reverse")
        # Conversion Tabular <==>Fasta
        DictInfo["FastaRenamedForward"] = MappingClass.tab2fasta(
            DictInfo["TabularRenamedForward"], "_forward")
        DictInfo["FastaRenamedReverse"] = MappingClass.tab2fasta(
            DictInfo["TabularRenamedReverse"], "_reverse")
        # Merging Pairs
        DictInfo["Trimmedreadconcatenated"] = MappingClass.concatenateforrev(
            [DictInfo["FastaRenamedForward"], DictInfo["FastaRenamedReverse"]])
        # Mapping steps for Single-End
        DictInfo["blastoutput"] = MappingClass.callmultiblastn(
            fasta=options.fastasequence,
            multifasta=DictInfo["Trimmedreadconcatenated"],
            outputformat=outformat7,
            suffix='_blastn.txt')
    elif (InpClass.sequencingtype == "Paired-End") and (InpClass.readforwardtype == "fasta"):
        if (InpClass.primer5forward is not None) and (InpClass.primer5reverse is not None):
            DictInfo["Trimmed5paired"] = TrimPaired.trimming5paired()
            # Conversion  Fasta<==>Tabular
            DictInfo["TabularReadsForward"] = MappingClass.fasta2tabular(
                imp=DictInfo["Trimmed5paired"][0], prefix="forward")
            DictInfo["TabularReadsReverse"] = MappingClass.fasta2tabular(
                imp=DictInfo["Trimmed5paired"][1], prefix="reverse")
        elif (InpClass.primer5forward is None) and (InpClass.primer3forward is None) and (InpClass.primer5reverse is None) and (InpClass.primer3reverse is None):
            # Conversion Fastq<==>Fasta Paired-end
            DictInfo["TabularReadsForward"] = MappingClass.fasta2tabular(
                imp=InpClass.readforward, prefix="forward")
            DictInfo["TabularReadsReverse"] = MappingClass.fasta2tabular(
                imp=InpClass.readreverse, prefix="reverse")
        else:
            log.write(msg116)
            sys.exit(1)
        # Rename Sequence ID forward
        DictInfo["TabularRenamedForward"] = MappingClass.seqrename(
            tabular=DictInfo["TabularReadsForward"], readirection="forward")
        DictInfo["TabularRenamedReverse"] = MappingClass.seqrename(
            tabular=DictInfo["TabularReadsReverse"], readirection="reverse")
        # Conversion Tabular <==>Fasta
        DictInfo["FastaRenamedForward"] = MappingClass.tab2fasta(
            DictInfo["TabularRenamedForward"], "_forward")
        DictInfo["FastaRenamedReverse"] = MappingClass.tab2fasta(
            DictInfo["TabularRenamedReverse"], "_reverse")
        # Merging Pairs
        DictInfo["Trimmedreadconcatenated"] = MappingClass.concatenateforrev(
            [DictInfo["FastaRenamedForward"], DictInfo["FastaRenamedReverse"]])
        # Mapping steps for Single-End
        DictInfo["blastoutput"] = MappingClass.callmultiblastn(
            fasta=options.fastasequence,
            multifasta=DictInfo["Trimmedreadconcatenated"],
            outputformat=outformat7,
            suffix='_blastn.txt')
    # End Paired fastq
    DictInfo["blastoutputnohash"] = MappingClass.hashclean(
        blastnout=DictInfo["blastoutput"], prefix="_blastn_nohash")
    # Filter reads steps (NO open-gaps, mismatch)
    DictInfo["blastoutputnohashfiltered"] = MappingClass.blastnfiltering(
        blastnout=DictInfo["blastoutputnohash"])
    MappingClass.cleantempfile()
