import optparse
from pyinteraseq_trimming import *
from pyinteraseq_mapping import *
from pytinteraseq_domains_definition import *
from pyinteraseq_genomefileparsing import *


parser = optparse.OptionParser(usage='python %prog Main PyInteract', version='1.0',)
parser.add_option('--readforward', action="store", dest="readforward", default=None,
                  help='Read dataset input forward')
parser.add_option('--readreverse', action="store", dest="readreverse", default=None,
                  help='Read dataset input reverse',)
parser.add_option('--readforwardtype', type='choice', choices=['fastq', 'fasta'], default=None,
                  help='Read dataset input forward')
parser.add_option('--readreversetype', type='choice', choices=['fastq', 'fasta'], default=None,
                  help='Read dataset input reverse')
parser.add_option('--sampletype', type='choice', choices=['genomic', 'selection', 'control'], default=None,
                  help='Select type of dataset')

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
reference_opts.add_option('--annotation', action="store", dest="annotation", default=None,
                          help='Annotation File(.gff|.bed)')
reference_opts.add_option('--chromosomename', action="store", dest="chromosomename", default=None,
                          help='Chromosome Name.(NC_XXXX)')
parser.add_option_group(reference_opts)

parser.add_option_group(query_opts)

reference_opts = optparse.OptionGroup(
    parser, 'Advanced Options',
    'Options for advanced analysis.',
    )
reference_opts.add_option('--minclonelength', action="store", dest="minclonelength", default='100',
                          help='Minumum clones length.')
reference_opts.add_option('--thread', action="store", dest="thread", default='1',
                          help='Number of thread.')
parser.add_option_group(reference_opts)

options, args = parser.parse_args()

if __name__ == '__main__':
    DictInfo = dict()
    DictFile = dict()
    outformat6 = '6 sseqid sstart send qseqid score sstrand'
    outformat7 = '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq'
    MergeFileList = []
    if options.readreverse is not None:
        if (options.readforwardtype == "fastq") and (options.readreversetype == "fastq"):
            DictInfo.update({
                "LogFilePath": TrimmingPaired(optparseinstance=options).logfilecreation(),
                "CutadaptPath": TrimmingSingle(optparseinstance=options).cutadaptchech(),
                "PickOtus": TrimmingSingle(optparseinstance=options).pickotuscheck(),
                "PickRepSeq": TrimmingSingle(optparseinstance=options).pickrepseqcheck(),
                "LogInfoAppended": TrimmingPaired(optparseinstance=options).inputinformationappen()
            })
    else:
        DictInfo.update({
            "LogFilePath": TrimmingSingle(optparseinstance=options).logfilecreation(),
            "CutadaptPath": TrimmingSingle(optparseinstance=options).cutadaptchech(),
            "PickOtus": TrimmingSingle(optparseinstance=options).pickotuscheck(),
            "PickRepSeq": TrimmingSingle(optparseinstance=options).pickrepseqcheck(),
            "LogInfoAppended": TrimmingSingle(optparseinstance=options).inputinformationappen()
        })
    DictInfo["fasta"] = GenomeFile(optparseinstance=options).fastareference()
    #
    DictInfo["annotation"] = AnnotationFile(optparseinstance=options).annotationbuild()
    #
    if options.readforwardtype == 'fastq':
        if options.readreverse is not None:
            DictInfo["Trimmed5paired"] = TrimmingPaired(optparseinstance=options).trimming5paired()
        else:
            DictInfo["Trimmed5single"] = TrimmingSingle(optparseinstance=options).trimming5single()

        if options.readreverse is not None:
            DictInfo["FastaReadsForward"] = BlastNlucleotide(optparseinstance=options).fastq2fasta(
                fastq=DictInfo["Trimmed5paired"][0], nameid="forward")
            #
            DictInfo["FastaReadsReverse"] = BlastNlucleotide(optparseinstance=options).fastq2fasta(
                fastq=DictInfo["Trimmed5paired"][1], nameid="reverse")
        else:
            DictInfo["FastaReadsForward"] = BlastNlucleotide(optparseinstance=options).fastq2fasta(
                fastq=DictInfo["Trimmed5single"], nameid="forward")
        DictInfo["TabularReadsForward"] = BlastNlucleotide(optparseinstance=options).fasta2tabular(
            imp=DictInfo["FastaReadsForward"], prefix="forward")
        if options.readreverse is not None:
            DictInfo["TabularReadsReverse"] = BlastNlucleotide(optparseinstance=options).fasta2tabular(
                imp=DictInfo["FastaReadsReverse"], prefix="reverse")
        DictInfo["TabularRenamedForward"] = BlastNlucleotide(optparseinstance=options).seqrename(
            DictInfo["TabularReadsForward"], "forward")
        if options.readreverse is not None:
            DictInfo["TabularRenamedReverse"] = BlastNlucleotide(optparseinstance=options).seqrename(
                DictInfo["TabularReadsReverse"], "reverse")
        DictInfo["FastaRenamedForward"] = BlastNlucleotide(optparseinstance=options).tab2fasta(
            DictInfo["TabularRenamedForward"], "_forward")

        if options.readreverse is not None:
            DictInfo["FastaRenamedReverse"] = BlastNlucleotide(optparseinstance=options).tab2fasta(
                DictInfo["TabularRenamedReverse"], "_reverse")
            DictInfo["Trimmedreadconcatenated"] = TrimmingPaired(optparseinstance=options).concatenateforrev(
                [DictInfo["FastaRenamedForward"], DictInfo["FastaRenamedReverse"]])
            DictInfo["blastoutput"] = BlastNlucleotide(optparseinstance=options).callmultiblastn(
                fasta=DictInfo["fasta"],
                multifasta=DictInfo["Trimmedreadconcatenated"],
                outputformat=outformat7,
                suffix='_blastn.txt')
        else:
            DictInfo["blastoutput"] = BlastNlucleotide(optparseinstance=options).callmultiblastn(
                fasta=DictInfo["fasta"],
                multifasta=DictInfo["FastaRenamedForward"],
                outputformat=outformat7,
                suffix='_blastn.txt')
        DictInfo["blastoutputnohash"] = DomainsDefinition(optparseinstance=options).hashclean(
            DictInfo["blastoutput"], "_blastn_nohash")
        DictInfo["blastoutputnohashfiltered"] = DomainsDefinition(optparseinstance=options).blastnfiltering(
            DictInfo["blastoutputnohash"])
        if options.readreverse is not None:
            DictInfo["blastoutputnohashfilteredfasta"] = BlastNlucleotide(optparseinstance=options).tab2fasta(
                DictInfo["blastoutputnohashfiltered"], "_paired")
            DictInfo["clustering"] = DomainsDefinition(optparseinstance=options).clustering(
                DictInfo["blastoutputnohashfilteredfasta"], "_paired")
        else:
            DictInfo["blastoutputnohashfilteredfasta"] = BlastNlucleotide(optparseinstance=options).tab2fasta(
                DictInfo["blastoutputnohashfiltered"], "_single")
            DictInfo["clustering"] = DomainsDefinition(optparseinstance=options).clustering(
                DictInfo["blastoutputnohashfilteredfasta"], "_single")
        #
        DictInfo["pickedreads"] = DomainsDefinition(optparseinstance=options).pickrepseq(
            DictInfo["clustering"], DictInfo["blastoutputnohashfilteredfasta"])
        DictInfo["pickedreadscleand"] = DomainsDefinition(optparseinstance=options).pysed(
            DictInfo["pickedreads"], '_clean.fasta', '-', '')
        #
        # DictInfo["blastedclones"] = BlastNlucleotide(optparseinstance=options).blastnclones(
        #     DictInfo["pickedreadscleand"], DictInfo["fasta"])
        DictInfo["blastedclones"] = BlastNlucleotide(optparseinstance=options).callmultiblastn(
            fasta=DictInfo["fasta"],
            multifasta=DictInfo["pickedreadscleand"],
            outputformat=outformat6,
            suffix='_blastnclones.tab')
        #
        DictInfo["bedparsed"] = DomainsDefinition(optparseinstance=options).bedparsing(DictInfo["blastedclones"])
        #
        DictInfo["clonesannotated"] = DomainsDefinition(optparseinstance=options).bedtoolsannotatefiltering(
            DomainsDefinition(optparseinstance=options).bedtoolsannotate(
                DictInfo["bedparsed"], DictInfo["annotation"]), 0.7)
        #
        DictInfo["clustercount"] = DomainsDefinition(optparseinstance=options).clonescount(
            DictInfo["clustering"])
        #
        DictInfo["clonescounted"] = DomainsDefinition(optparseinstance=options).mergingcount(
            DictInfo["bedparsed"], DictInfo["clustercount"])
        #
        #
        DictInfo["clonescountedfiltered"] = DomainsDefinition(optparseinstance=options).filteringclonescount(
            DictInfo["clonescounted"], 10)
        #
        DictInfo["clonescountedmerged"] = DomainsDefinition(optparseinstance=options).pybedtoolsmerge(
            DictInfo["clonescountedfiltered"])
        #
        DictInfo["clonesmergedfasta"] = DomainsDefinition(optparseinstance=options).pybedtoolstofasta(
            DictInfo["clonescountedmerged"], DictInfo["fasta"])
        #
        DictInfo["tabwithdescription"] = DomainsDefinition(optparseinstance=options).adddescription(
            DictInfo["clonescountedmerged"], DictInfo["annotation"], 0.7)
        DictInfo["clonenseqfasta"] = BlastNlucleotide(optparseinstance=options).fasta2tabular(
            imp=DictInfo["clonesmergedfasta"], prefix='_clonestabular')
        DictFile["tabwithsequence"] = DomainsDefinition(optparseinstance=options).addsequence(
            outputfromdescription=DictInfo["tabwithdescription"], outputfasta2tab=DictInfo["clonenseqfasta"])
        print DictInfo
        DomainsDefinition(optparseinstance=options).cleantemporaryfilesinglend(DictInfo)
    elif options.readforwardtype == 'fasta':
        if options.readreverse is not None:
            DictInfo["Trimmed5paired"] = TrimmingPaired(optparseinstance=options).trimming5paired()
        else:
            DictInfo["Trimmed5single"] = TrimmingSingle(optparseinstance=options).trimming5single()

        if options.readreverse is not None:
            DictInfo["FastaReadsForward"] = BlastNlucleotide(optparseinstance=options).fasta2tabular(
                imp=DictInfo["Trimmed5paired"][0], prefix="forward")
            #
            DictInfo["FastaReadsReverse"] = BlastNlucleotide(optparseinstance=options).fasta2tabular(
                imp=DictInfo["Trimmed5paired"][1], prefix="reverse")
        else:
            DictInfo["FastaReadsForward"] = BlastNlucleotide(optparseinstance=options).fasta2tabular(
                imp=DictInfo["Trimmed5single"], prefix="forward")
        DictInfo["TabularRenamedForward"] = BlastNlucleotide(optparseinstance=options).seqrename(
            DictInfo["FastaReadsForward"], "forward")
        if options.readreverse is not None:
            DictInfo["TabularRenamedReverse"] = BlastNlucleotide(optparseinstance=options).seqrename(
                DictInfo["FastaReadsReverse"], "reverse")
        DictInfo["FastaRenamedForward"] = BlastNlucleotide(optparseinstance=options).tab2fasta(
            DictInfo["TabularRenamedForward"], "_forward")
        if options.readreverse is not None:
            DictInfo["FastaRenamedReverse"] = BlastNlucleotide(optparseinstance=options).tab2fasta(
                DictInfo["TabularRenamedReverse"], "_reverse")
            DictInfo["Trimmedreadconcatenated"] = TrimmingPaired(optparseinstance=options).concatenateforrev(
                [DictInfo["FastaRenamedForward"], DictInfo["FastaRenamedReverse"]])
            DictInfo["blastoutput"] = BlastNlucleotide(optparseinstance=options).callmultiblastn(
                fasta=DictInfo["fasta"],
                multifasta=DictInfo["Trimmedreadconcatenated"],
                outputformat=outformat7,
                suffix='_blastn.txt')
        else:
            DictInfo["blastoutput"] = BlastNlucleotide(optparseinstance=options).callmultiblastn(
                fasta=DictInfo["fasta"],
                multifasta=DictInfo["FastaRenamedForward"],
                outputformat=outformat7,
                suffix='_blastn.txt')
        DictInfo["blastoutputnohash"] = DomainsDefinition(optparseinstance=options).hashclean(
            DictInfo["blastoutput"], "_blastn_nohash")
        DictInfo["blastoutputnohashfiltered"] = DomainsDefinition(optparseinstance=options).blastnfiltering(
            DictInfo["blastoutputnohash"])
        if options.readreverse is not None:
            DictInfo["blastoutputnohashfilteredfasta"] = BlastNlucleotide(optparseinstance=options).tab2fasta(
                DictInfo["blastoutputnohashfiltered"], "_paired")
            DictInfo["clustering"] = DomainsDefinition(optparseinstance=options).clustering(
                DictInfo["blastoutputnohashfilteredfasta"], "_paired")
        else:
            DictInfo["blastoutputnohashfilteredfasta"] = BlastNlucleotide(optparseinstance=options).tab2fasta(
                DictInfo["blastoutputnohashfiltered"], "_single")
            DictInfo["clustering"] = DomainsDefinition(optparseinstance=options).clustering(
                DictInfo["blastoutputnohashfilteredfasta"], "_single")
        #
        DictInfo["pickedreads"] = DomainsDefinition(optparseinstance=options).pickrepseq(
            DictInfo["clustering"], DictInfo["blastoutputnohashfilteredfasta"])
        DictInfo["pickedreadscleand"] = DomainsDefinition(optparseinstance=options).pysed(
            DictInfo["pickedreads"], '_clean.fasta', '-', '')
        #
        DictInfo["blastedclones"] = BlastNlucleotide(optparseinstance=options).callmultiblastn(
            fasta=DictInfo["fasta"],
            multifasta=DictInfo["pickedreadscleand"],
            outputformat=outformat6,
            suffix='_blastnclones.tab')
        #
        DictInfo["bedparsed"] = DomainsDefinition(optparseinstance=options).bedparsing(DictInfo["blastedclones"])
        #
        DictInfo["clonesannotated"] = DomainsDefinition(optparseinstance=options).bedtoolsannotatefiltering(
            DomainsDefinition(optparseinstance=options).bedtoolsannotate(
                DictInfo["bedparsed"], DictInfo["annotation"]), 0.7)
        #
        DictInfo["clustercount"] = DomainsDefinition(optparseinstance=options).clonescount(
            DictInfo["clustering"])
        # #
        DictInfo["clonescounted"] = DomainsDefinition(optparseinstance=options).mergingcount(
            DictInfo["bedparsed"], DictInfo["clustercount"])
        #
        DictInfo["clonescountedfiltered"] = DomainsDefinition(optparseinstance=options).filteringclonescount(
            DictInfo["clonescounted"], 10)
        #
        DictInfo["clonescountedmerged"] = DomainsDefinition(optparseinstance=options).pybedtoolsmerge(
            DictInfo["clonescountedfiltered"])
        #
        DictInfo["clonesmergedfasta"] = DomainsDefinition(optparseinstance=options).pybedtoolstofasta(
            DictInfo["clonescountedmerged"], DictInfo["fasta"])
        #
        DictInfo["tabwithdescription"] = DomainsDefinition(optparseinstance=options).adddescription(
            DictInfo["clonescountedmerged"], DictInfo["annotation"], 0.7)
        DictInfo["clonenseqfasta"] = BlastNlucleotide(optparseinstance=options).fasta2tabular(
            imp=DictInfo["clonesmergedfasta"], prefix='_clonestabular')
        DictFile["tabwithsequence"] = DomainsDefinition(optparseinstance=options).addsequence(
            outputfromdescription=DictInfo["tabwithdescription"], outputfasta2tab=DictInfo["clonenseqfasta"])
        print DictInfo
        DomainsDefinition(optparseinstance=options).cleantemporaryfilesinglend(DictInfo)

