import optparse
from pyinteraseq_mapping import *
from output_message import *
import sys
import traceback

parser = optparse.OptionParser(usage='python %prog pyinteraseq_main_mapping.py', version='1.0',)
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
                       help='Primer 5\'forward')
primer_opts.add_option('--primer3forward', action="store", dest="primer3forward", default=None,
                       help='Primer 3\'forward')
primer_opts.add_option('--primer5reverse', action="store", dest="primer5reverse", default=None,
                       help='Primer 5\'reverse')
primer_opts.add_option('--primer3reverse', action="store", dest="primer3reverse", default=None,
                       help='Primer 3\'reverse')
parser.add_option_group(primer_opts)

reference_opts = optparse.OptionGroup(
    parser, 'Primers Options',
    'Options for primer sequence and trimming.',
    )
reference_opts.add_option('--fastasequence', action="store", dest="fastasequence", default=None,
                          help='Genome sequence fasta file.(.fasta|.fna|.fa)')
reference_opts.add_option('--annotationpath', action="store", dest="annotationpath", default=None,
                          help='annotation path')
reference_opts.add_option('--thread', action="store", dest="thread", default='10',
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
advance_opts.add_option('--mismatch', action="store", dest="mismatch", type="float",
                        default=3.0, help='Percentage of Mismatch allowed.')
parser.add_option_group(advance_opts)

options, args = parser.parse_args()

if __name__ == '__main__':
    DictInfo = dict()
    MappingClass = BlastNlucleotide(optparseinstance=options)
    if MappingClass.sequencingtype == "Paired-End":
        if (MappingClass.readforwardtype is "fastq") and (MappingClass.readreversetype is "fastq"):
            DictInfo.update({
                "LogInfoAppended": MappingClass.inputinformationappen()
            })
        elif (MappingClass.readforwardtype is "fasta") and (MappingClass.readreversetype is "fasta"):
            DictInfo.update({
                "LogInfoAppended": MappingClass.inputinformationappen()
            })
        else:
            log = open(MappingClass.inputfilelog, "a", 0)
            log.write(msg5)
            sys.exit(1)
    elif MappingClass.sequencingtype == "Single-End":
        DictInfo.update({
            "LogInfoAppended": MappingClass.inputinformationappen()
        })
    log = open(MappingClass.inputfilelog, "a", 0)
    if MappingClass.sequencingtype == "Single-End":
        log.write(msg127)
        try:
            if (MappingClass.primer5forward is not None) and (MappingClass.primer3forward is not None):
                DictInfo["forward5trimmed"] = MappingClass.trimming5single(readfile=MappingClass.readforward,
                                                                           primer5=MappingClass.primer5forward,
                                                                           direction="_forward5trimmed.",
                                                                           readtype=MappingClass.readforwardtype)
                DictInfo["Trimmedreadconcatenated"] = MappingClass.trimming3single(readfile=DictInfo["forward5trimmed"],
                                                                                   primer3=MappingClass.primer3forward,
                                                                                   direction="_forward3trimmed.",
                                                                                   readtype=MappingClass.readforwardtype)
            elif (MappingClass.primer5forward is None) and (MappingClass.primer3forward is None):
                DictInfo["Trimmedreadconcatenated"] = MappingClass.readforward
            else:
                DictInfo["Trimmedreadconcatenated"] = MappingClass.trimming5single(readfile=MappingClass.readforward,
                                                                                   primer5=MappingClass.primer5forward,
                                                                                   direction="_forward5trimmed.",
                                                                                   readtype=MappingClass.readforwardtype)
        except traceback:
            log.write(msg41)
            sys.exit(1)
        else:
            log.write(msg130)
    elif MappingClass.sequencingtype == "Paired-End":
        log.write(msg129)
        try:
            if (MappingClass.primer5forward is not None) and (MappingClass.primer3forward is not None) and (
                        MappingClass.primer5reverse is not None) and (MappingClass.primer3reverse is not None):
                DictInfo["forward5trimmed"] = MappingClass.trimming5single(readfile=MappingClass.readforward,
                                                                           primer5=MappingClass.primer5forward,
                                                                           direction="_forward5trimmed.",
                                                                           readtype=MappingClass.readforwardtype)
                DictInfo["reverse5trimmed"] = MappingClass.trimming5single(readfile=MappingClass.readreverse,
                                                                           primer5=MappingClass.primer5reverse,
                                                                           direction="_reverse5trimmed.",
                                                                           readtype=MappingClass.readreversetype)
                DictInfo["forward3trimmed"] = MappingClass.trimming3single(readfile=DictInfo["forward5trimmed"],
                                                                           primer3=MappingClass.primer3forward,
                                                                           direction="_forward3trimmed.",
                                                                           readtype=MappingClass.readforwardtype)
                DictInfo["reverse3trimmed"] = MappingClass.trimming3single(readfile=DictInfo["reverse5trimmed"],
                                                                           primer3=MappingClass.primer3reverse,
                                                                           direction="_reverse3trimmed.",
                                                                           readtype=MappingClass.readforwardtype)
                DictInfo["Trimmedreadconcatenated"] = MappingClass.concatenateforrev([DictInfo["forward3trimmed"],
                                                                                      DictInfo["reverse3trimmed"]])
            elif (MappingClass.primer5forward is None) and (MappingClass.primer3forward is None) and (
                        MappingClass.primer5reverse is None) and (MappingClass.primer3reverse is None):
                DictInfo["Trimmedreadconcatenated"] = MappingClass.concatenateforrev([MappingClass.readforward,
                                                                                      MappingClass.readreverse])
            elif (MappingClass.primer5forward is not None) and (MappingClass.primer3forward is not None) and (
                        MappingClass.primer5reverse is None) and (MappingClass.primer3reverse is None):
                DictInfo["forward5trimmed"] = MappingClass.trimming5single(readfile=MappingClass.readforward,
                                                                           primer5=MappingClass.primer5forward,
                                                                           direction="_forward5trimmed.",
                                                                           readtype=MappingClass.readforwardtype)
                DictInfo["forward3trimmed"] = MappingClass.trimming3single(readfile=DictInfo["forward5trimmed"],
                                                                           primer3=MappingClass.primer3forward,
                                                                           direction="_forward3trimmed.",
                                                                           readtype=MappingClass.readforwardtype)
                DictInfo["Trimmedreadconcatenated"] = MappingClass.concatenateforrev([DictInfo["forward3trimmed"],
                                                                                      MappingClass.readreverse])
            elif (MappingClass.primer5forward is None) and (MappingClass.primer3forward is None) and (
                        MappingClass.primer5reverse is not None) and (MappingClass.primer3reverse is not None):
                DictInfo["reverse5trimmed"] = MappingClass.trimming5single(readfile=MappingClass.readreverse,
                                                                           primer5=MappingClass.primer5reverse,
                                                                           direction="_reverse5trimmed.",
                                                                           readtype=MappingClass.readreversetype)
                DictInfo["reverse3trimmed"] = MappingClass.trimming3single(readfile=DictInfo["reverse5trimmed"],
                                                                           primer3=MappingClass.primer3reverse,
                                                                           direction="_reverse3trimmed.",
                                                                           readtype=MappingClass.readforwardtype)
            elif (MappingClass.primer5forward is not None) and (MappingClass.primer3forward is None) and (
                        MappingClass.primer5reverse is not None) and (MappingClass.primer3reverse is None):
                DictInfo["reverse5trimmed"] = MappingClass.trimming5single(readfile=MappingClass.readreverse,
                                                                           primer5=MappingClass.primer5reverse,
                                                                           direction="_reverse5trimmed.",
                                                                           readtype=MappingClass.readreversetype)
                DictInfo["forward5trimmed"] = MappingClass.trimming5single(readfile=MappingClass.readforward,
                                                                           primer5=MappingClass.primer5forward,
                                                                           direction="_forward5trimmed.",
                                                                           readtype=MappingClass.readforwardtype)                           
            else:
                log.write(msg4)
                sys.exit(1)
        except traceback:
            log.write(msg41)
            sys.exit(1)
	else:
            log.write(msg130)

    #MappingClass.gtf = options.annotationpath + ".chr.gtf"
    MappingClass.gtf = "/home/spuccio/PhageRnaBinding/Homo_sapiens.GRCh38.99.chr.gtf"
    MappingClass.genomesize = "/home/spuccio/PhageRnaBinding/sizes.genome"
    MappingClass.genomeorganism = "/home/spuccio/PhageRnaBinding/HStranscriptome.idx"

    if MappingClass.sequencingtype == "Single-End":
    	DictInfo["indexvalues"] = MappingClass.kallistoindexvalues(DictInfo["Trimmedreadconcatenated"])
	DictInfo["bam"] = MappingClass.mappingkallisto(MappingClass.genomeorganism,MappingClass.gtf,MappingClass.genomesize,DictInfo["indexvalues"][0],DictInfo["indexvalues"][1],DictInfo["Trimmedreadconcatenated"])
	DictInfo["samtr"] = MappingClass.mappingkallistotran(MappingClass.genomeorganism, DictInfo["indexvalues"][0],DictInfo["indexvalues"][1],DictInfo["Trimmedreadconcatenated"])
	DictInfo["sortedbam"] = MappingClass.sortbam(bamfile=DictInfo["samtr"])
    elif MappingClass.sequencingtype == "Paired-End":
        DictInfo["bam"] = MappingClass.mappingkallistopaired(MappingClass.genomeorganism,MappingClass.gtf,MappingClass.genomesize,DictInfo["forward5trimmed"],DictInfo["reverse5trimmed"])
	DictInfo["samtr"] = MappingClass.mappingkallistopairedttran(MappingClass.genomeorganism,DictInfo["forward5trimmed"],DictInfo["reverse5trimmed"])
    	DictInfo["sortedbam"] = MappingClass.sortbam(bamfile=DictInfo["samtr"])

    DictInfo["sam"] = MappingClass.conversionbam2sam(bamfile=DictInfo["bam"])
    DictInfo["filtsam"] = MappingClass.filtersam(samfile=DictInfo["sam"])
    DictInfo["filtbam"] = MappingClass.conversionsam2bam(samfile=DictInfo["filtsam"])
    DictInfo["bedgraph"] = MappingClass.bam2bedgraph(bamfile=DictInfo["filtbam"])
    DictInfo["bw"] = MappingClass.bedgraph2bw(bedgraph=DictInfo["bedgraph"])
    MappingClass.cleantempfile()
