import optparse
from pyinteraseq_trimming import *
from pyinteraseq_mapping import *
from pytinteraseq_domains_definition import *
from pyinteraseq_genomefileparsing import *
from Bio import SeqIO
from Bio.Alphabet import IUPAC

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
                          help='Genome sequence fasta file.')
reference_opts.add_option('--annotation', action="store", dest="annotation", default=None,
                          help='Annotation File')
parser.add_option_group(reference_opts)

parser.add_option_group(query_opts)

reference_opts = optparse.OptionGroup(
    parser, 'Advanced Options',
    'Options for advanced analysis.',
    )
reference_opts.add_option('--minclonelength', action="store", dest="minclonelength", default='100',
                          help='Minumum clones length.')
reference_opts.add_option('--thread', action="store", dest="thread", default='2',
                          help='Number of thread.')
parser.add_option_group(reference_opts)

options, args = parser.parse_args()

if __name__ == '__main__':
    if options.readreverse is not None:
        if (options.readforwardtype == "fastq") and (options.readreversetype == "fastq"):
            DictInfo = dict()
            MergeFileList = []
            DictInfo.update({
                "LogFilePath": TrimmingPaired(optparseinstance=options).logfilecreation(),
                "CutadaptPath": TrimmingSingle(optparseinstance=options).cutadaptchech(),
                "PickOtus": TrimmingSingle(optparseinstance=options).pickotuscheck(),
                "PickRepSeq": TrimmingSingle(optparseinstance=options).pickrepseqcheck(),
                "LogInfoAppended": TrimmingPaired(optparseinstance=options).inputinformationappen()
            })
            DictInfo["Trimmed5paired"] = TrimmingPaired(optparseinstance=options).trimming5paired()
            # DictInfo["Trimmed3paired"] = TrimmingPaired(optparseinstance=options).
            # trimming3paired()
            # DictInfo["MergedRead"] = TrimmingPaired(optparseinstance=options).
            # concatenateforrev(readlist=DictInfo["Trimmed3paired"])
            DictInfo["FastaReadsForward"] = BlastNlucleotide(optparseinstance=options).\
                fastq2fasta(fastq=DictInfo["Trimmed5paired"][0], nameid="forward")
            DictInfo["FastaReadsReverse"] = BlastNlucleotide(optparseinstance=options).\
                fastq2fasta(fastq=DictInfo["Trimmed5paired"][1], nameid="reverse")
            DictInfo["TabularReadsForward"] = BlastNlucleotide(optparseinstance=options).\
                fasta2tab(fasta=DictInfo["FastaReadsForward"], nameid="forward")
            DictInfo["TabularReadsReverse"] = BlastNlucleotide(optparseinstance=options).\
                fasta2tab(fasta=DictInfo["FastaReadsReverse"], nameid="forward")
            DictInfo["TabularRenamedForward"] = BlastNlucleotide(optparseinstance=options).\
                seqrename(DictInfo["TabularReadsForward"], "forward")
            DictInfo["TabularRenamedReverse"] = BlastNlucleotide(optparseinstance=options).\
                seqrename(DictInfo["TabularReadsReverse"], "reverse")
            DictInfo["FastaRenamedForward"] = BlastNlucleotide(optparseinstance=options).\
                tab2fasta(DictInfo["TabularRenamedForward"], "forward")
            DictInfo["FastaRenamedReverse"] = BlastNlucleotide(optparseinstance=options).\
                tab2fasta(DictInfo["TabularRenamedReverse"], "reverse")
            MergeFileList.append(DictInfo["FastaRenamedForward"])
            MergeFileList.append(DictInfo["FastaRenamedReverse"])
            # MergeFileList.append('/Users/simone/output_test/26695_1_newid.fasta')
            # MergeFileList.append('/Users/simone/output_test/26695_2_newid.fasta')
            #DictInfo["FastaMerged"] = BlastNlucleotide(optparseinstance=options).concatenateforrev(MergeFileList)
            #DictInfo["blastoutput"] = BlastNlucleotide(optparseinstance=options).blastn(DictInfo["FastaMerged"])
            #DictInfo["blastoutputnohash"] = DomainsDefinition(optparseinstance=options).hashclean(DictInfo["blastoutput"])

    else:
        DictInfo = dict()
        MergeFileList = []
        DictInfo.update({
            "LogFilePath": TrimmingSingle(optparseinstance=options).logfilecreation(),
            "CutadaptPath": TrimmingSingle(optparseinstance=options).cutadaptchech(),
            "PickOtus": TrimmingSingle(optparseinstance=options).pickotuscheck(),
            "PickRepSeq": TrimmingSingle(optparseinstance=options).pickrepseqcheck(),
            "LogInfoAppended": TrimmingSingle(optparseinstance=options).inputinformationappen()
        })
        if options.readforwardtype == "fastq":
            DictInfo["fasta"] = GenomeFile(optparseinstance=options).fastareference()
            DictInfo["annotation"] = AnnotationFile(optparseinstance=options).annotationbuild()
            # DictInfo["Trimmed5single"] = TrimmingSingle(optparseinstance=options).trimming5single()
            # DictInfo["FastaReadsForward"] = BlastNlucleotide(optparseinstance=options).\
            #     fastq2fasta(fastq=DictInfo["Trimmed5single"], nameid="forward")
            # DictInfo["TabularReadsForward"] = BlastNlucleotide(optparseinstance=options).\
            #     fasta2tab(fasta=DictInfo["FastaReadsForward"], nameid="forward")
            # DictInfo["TabularRenamedForward"] = BlastNlucleotide(optparseinstance=options).\
            #     seqrename(DictInfo["TabularReadsForward"], "forward")
            # DictInfo["FastaRenamedForward"] = BlastNlucleotide(optparseinstance=options).\
            #     tab2fasta(DictInfo["TabularRenamedForward"], "forward")
            # DictInfo["blastoutput"] = BlastNlucleotide(optparseinstance=options).blastn(DictInfo["FastaRenamedForward"],
            #                                                                             DictInfo["fasta"])
            # DictInfo["blastoutputnohash"] = DomainsDefinition(optparseinstance=options).hashclean(DictInfo["blastoutput"])
            # DictInfo["blastoutputnohashfiltered"] = DomainsDefinition(optparseinstance=options).blastnfiltering(DictInfo["blastoutputnohash"])
            # DictInfo["blastoutputnohashfilteredfasta"] = BlastNlucleotide(optparseinstance=options).tab2fasta(DictInfo["blastoutputnohashfiltered"], "single")
            # DictInfo["clustering"] = DomainsDefinition(optparseinstance=options).clustering(DictInfo["blastoutputnohashfilteredfasta"])
            # DictInfo["pickedreads"] = DomainsDefinition(optparseinstance=options).pickrepseq(DictInfo["clustering"],
            #                                                                                  DictInfo["blastoutputnohashfilteredfasta"])
            # DictInfo["pickedreadscleand"] = DomainsDefinition(optparseinstance=options).pysed(DictInfo["pickedreads"],
            #                                                                                   '_clean.fasta', '-', '')
            # DictInfo["blastedclones"] = BlastNlucleotide(optparseinstance=options).blastnclones(
            #     DictInfo["pickedreadscleand"], DictInfo["fasta"])
            # DictInfo["clustering"] = '/Users/simone/output_test/26695_picked/26695_newid_otus.txt'
            # DictInfo["blastedclones"] = '/Users/simone/output_test/26695_blastnclones.tab'
            # DictInfo["bedparsed"] = DomainsDefinition(optparseinstance=options).bedparsing(DictInfo["blastedclones"])
            # DictInfo["clonesannotated"] = DomainsDefinition(optparseinstance=options).bedtoolsannotatefiltering(
            #     DomainsDefinition(optparseinstance=options).bedtoolsannotate(
            #         DictInfo["bedparsed"], DictInfo["annotation"]), 0.7)
            # DictInfo["clustercount"] = DomainsDefinition(optparseinstance=options).clonescount(DictInfo["clustering"])
            # DictInfo["clonescounted"] = DomainsDefinition(optparseinstance=options).mergingcount(
            #     DictInfo["bedparsed"], DictInfo["clustercount"])
            # # #go to domain
            # #
            # DictInfo["clonescountedfiltered"] = DomainsDefinition(optparseinstance=options).filteringclonescount(
            #     DictInfo["clonescounted"], 10)
            # #
            # DictInfo["clonescountedmerged"] = DomainsDefinition(optparseinstance=options).pybedtoolsmerge(
            #     DictInfo["clonescountedfiltered"])
            # #
            # DictInfo["clonesmergedfasta"] = DomainsDefinition(optparseinstance=options).pybedtoolstofasta(
            #     DictInfo["clonescountedmerged"], DictInfo["fasta"])
            # #
            # DictInfo["tabwithdescription"] = DomainsDefinition(optparseinstance=options).adddescription(
            #     DictInfo["clonescountedmerged"], DictInfo["annotation"], 0.7)
            #
            # DictInfo["clonenseqfasta"] = BlastNlucleotide(optparseinstance=options).fasta2tabular(
            #     imp=DictInfo["clonesmergedfasta"], prefix='_clonestabular')
            # DictInfo["tabwithsequence"] = DomainsDefinition(optparseinstance=options).addsequence(
            #     outputfromdescription=DictInfo["tabwithdescription"], outputfasta2tab=DictInfo["clonenseqfasta"])
            DictInfo["clustering"] = '/Users/simone/output_test/26695_picked/26695_newid_otus.txt'
            DictInfo["blastedclones"] = '/Users/simone/output_test/26695_blastnclones.tab'
            DictInfo["bedparsed"] = DomainsDefinition(optparseinstance=options).bedparsing(DictInfo["blastedclones"])
            DictInfo["clonesannotated"] = DomainsDefinition(optparseinstance=options).bedtoolsannotatefiltering(
                DomainsDefinition(optparseinstance=options).bedtoolsannotate(
                    DictInfo["bedparsed"], DictInfo["annotation"]), 0.7)
            DictInfo["clustercount"] = DomainsDefinition(optparseinstance=options).clonescount(DictInfo["clustering"])
            DictInfo["clonescounted"] = DomainsDefinition(optparseinstance=options).mergingcount(
                DictInfo["bedparsed"], DictInfo["clustercount"])
            # #go to domain
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
            DictInfo["tabwithsequence"] = DomainsDefinition(optparseinstance=options).addsequence(
                outputfromdescription=DictInfo["tabwithdescription"], outputfasta2tab=DictInfo["clonenseqfasta"])
            #inputfile = '/Users/simone/output_test/26695_blastnclonesmerge.fasta'
            for seq_record in SeqIO.parse(DictInfo["clonesmergedfasta"], "fasta", alphabet=IUPAC.ambiguous_dna):
                DictInfo["allframes"] = DomainsDefinition(optparseinstance=options).translatednaframes(
                    seq_record, DictInfo["clonesmergedfasta"])
            # print DictInfo
        elif options.readforwardtype == "fasta":
            print 'Single-End e fasta'






