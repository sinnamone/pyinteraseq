import optparse
from pyinteraseq_mapping import *
from pytinteraseq_domains_definition import *

parser = optparse.OptionParser(usage='python %prog Main PyInteract', version='1.0',)
parser.add_option('--readforwardtrimmed', action="store", dest="readforwardtrimmed", default=None,
                  help='Read dataset input forward')
parser.add_option('--readreversetrimmed', action="store", dest="readreversetrimmed", default=None,
                  help='Read dataset input reverse',)
parser.add_option('--readforwardtrimmedtype', type='choice', choices=['fastq', 'fasta'], default=None,
                  help='Read dataset input forward')
parser.add_option('--readreversetrimmedtype', type='choice', choices=['fastq', 'fasta'], default=None,
                  help='Read dataset input reverse')
parser.add_option('--sampletype', type='choice', choices=['genomic', 'selection', 'control'], default=None,
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
reference_opts.add_option('--overlapintersect', action="store", dest="overlapintersect", type="float",
                          default=0.7, help='Parameters -f of bedtools intersect.')
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
    DictInfo["blastoutputnohash"] = DomainsDefinition(optparseinstance=options).hashclean(
        blastnout=DictInfo["blastoutput"], prefix="_blastn_nohash")
    # Filter reads steps (NO open-gaps, mismatch)
    DictInfo["blastoutputnohashfiltered"] = DomainsDefinition(optparseinstance=options).blastnfiltering(
        blastnout=DictInfo["blastoutputnohash"])
    # Conversion Filter Blastn Table in Fasta
    DictInfo["blastoutputnohashfilteredfasta"] = BlastNlucleotide(optparseinstance=options).tab2fasta(
        tabular=DictInfo["blastoutputnohashfiltered"], prefixoutput="_blastnfiltered")
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
    DictInfo["blastedclones"] = BlastNlucleotide(optparseinstance=options).callmultiblastn(
        fasta=DictInfo["fasta"],
        multifasta=DictInfo["pickedreadscleand"],
        outputformat=outformat6,
        suffix='_blastnclones.tab')
    # Parser that take as input blastn format 6 and create a standard bed6
    DictInfo["bedparsed"] = DomainsDefinition(optparseinstance=options).bedparsing(DictInfo["blastedclones"])
    # Filtering the output of Bedtools annotate (call inside the function) using a flot percentage of overlap
    DictInfo["clonesannotated"] = DomainsDefinition(optparseinstance=options).bedtoolsannotatefiltering(
        # Bedtools annotate to identify the clones inside the CDS
        DomainsDefinition(optparseinstance=options).bedtoolsannotate(
            DictInfo["bedparsed"], DictInfo["annotation"]), options.overlapintersect)
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
        pybedtoolsmergeoutput=DictInfo["clonescountedmerged"], fastqsequence=DictInfo["fasta"])
    # Add description present in ptt 0.7 is the overlap of intersect cds
    DictInfo["tabwithdescription"] = DomainsDefinition(optparseinstance=options).adddescription(
        clonesmerged=DictInfo["clonescountedmerged"], annotation=DictInfo["annotation"],
        percthr=options.overlapintersect)
    # Conversion Fasta ==>Tabular
    DictInfo["clonenseqfasta"] = BlastNlucleotide(optparseinstance=options).fasta2tabular(
        imp=DictInfo["clonesmergedfasta"], prefix='_clonestabular')
    # Add sequence to output table
    DictFile["tabwithsequence"] = DomainsDefinition(optparseinstance=options).addsequence(
        outputfromdescription=DictInfo["tabwithdescription"], outputfasta2tab=DictInfo["clonenseqfasta"])
    # print DictInfo
    # DomainsDefinition(optparseinstance=options).cleantemporaryfilesinglend(DictInfo)
