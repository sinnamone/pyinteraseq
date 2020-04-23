import subprocess
from output_message import *
import os, sys
from pyinteraseq_trimming import Trimming
from multiprocessing import Pool
import traceback, urllib, collections, pysam
from Bio import SeqIO
import numpy as np

class BlastNlucleotide(Trimming):

    def __init__(self, optparseinstance):
        Trimming.__init__(self, optparseinstance)
        self.header = []
        self.pool = Pool(processes=int(self.thread))
        self.dbname = self.outputfolder + os.path.basename(self.fastasequence.split('/')[(-1)]).split('.')[0]
        self.filelog = None
        self.parts = ''
        #self.genomesize = "/Users/simone/Pyinteraseq/PyinteraseqEucaryotic/output/HS/sizes.genome"
        #self.kallistoindex="/Users/simone/Pyinteraseq/PyinteraseqEucaryotic/output/HS/Homo_sapiens"
        #self.gtf = "/Users/simone/Pyinteraseq/PyinteraseqEucaryotic/output/HS/Homo_sapiens.GRCh38.99.chr.gtf"
        #self.genomesize = "/home/spuccio/PhageRnaBinding/size.genome"
        #self.kallistoindex="/home/spuccio/PhageRnaBinding/HStranscriptome.idx"
        #self.gtf = "/home/spuccio/PhageRnaBinding//Homo_sapiens.GRCh38.99.chr.gtf"
        return


    def parsesamkallistofiles(self, samfile):
        """

        :param samfile:
        :return:
        """
        self.filelog = self.logopen()
        try:
            with open(samfile) as (infile):
                for line in infile:
                    if line.startswith('@'):
                        continue
                    if line.split('\t')[5] is '*':
                        continue
                    self.parts = line.split('\t')
                    assert len(self.parts) == len(self.samInfoFallistoFields)
                    normalizedinfo = {'QNAME': None if self.parts[0] == '.' else urllib.unquote(self.parts[0]), 
                       'FLAG': None if self.parts[1] == '.' else urllib.unquote(self.parts[1]), 
                       'RNAME': None if self.parts[2] == '.' else urllib.unquote(self.parts[2]), 
                       'POS': None if self.parts[3] == '.' else urllib.unquote(self.parts[3]), 
                       'MAPQ': None if self.parts[4] == '.' else urllib.unquote(self.parts[4]), 
                       'CIGAR': None if self.parts[5] == '.' else urllib.unquote(self.parts[5]), 
                       'RNEXT': None if self.parts[6] == '.' else urllib.unquote(self.parts[6]), 
                       'PNEXT': None if self.parts[7] == '.' else urllib.unquote(self.parts[7]), 
                       'TLEN': None if self.parts[8] == '.' else urllib.unquote(self.parts[8]), 
                       'SEQ': None if self.parts[9] == '.' else urllib.unquote(self.parts[9]), 
                       'QUAL': None if self.parts[10] == '.' else urllib.unquote(self.parts[10]), 
                       'NH': None if self.parts[10] == '.' else urllib.unquote(self.parts[11])}
                    yield self.SAMRecordKallisto(**normalizedinfo)

        except traceback:
            self.filelog.write(msg63)
            sys.exit(1)
        else:
            self.filelog.write(msg64)

        return


    def conversionbam2sam(self, bamfile):
        self.filelog = self.logopen()
        fnull = open(os.devnull, 'w')
        try:
            subprocess.check_call(['samtools', 'view', '-h', bamfile, '-o', ('/').join([self.out, ('').join([self.outputid, '.sam'])])], stderr=fnull, stdout=fnull)
        except traceback:
            self.filelog.write(msg82)
            sys.exit(1)
        else:
            self.filelog.write(msg83)
            return ('/').join([self.out, ('').join([self.outputid, '.sam'])])


    def filtersam(self, samfile):
        self.filelog = self.logopen()
        fnull = open(os.devnull, 'w')
        try:
            with open( ("/").join([self.out, ('').join([self.outputid, '_filtered.sam'])]), 'w') as filtsam:
            	subprocess.check_call(['nawk','''BEGIN{FS=OFS="\t"}{if($3!=\"*\") print $0} ''', samfile], stderr=fnull, stdout=filtsam)
        except traceback:
            self.filelog.write(msg82)
            sys.exit(1)
        else:
            self.filelog.write(msg83)
            filtsam.close()
            return ("/").join([self.out, ('').join([self.outputid, '_filtered.sam'])])

    def conversionsam2bam(self, samfile):
        self.filelog = self.logopen()
        try:
            pysam.view("-hbS", samfile, "-@", self.thread, "-o", ('/').join([self.out, ('').join([self.outputid, '_filtered.bam'])]), catch_stdout=False)
        except traceback:
            self.filelog.write(msg82)
            sys.exit(1)
        else:
            self.filelog.write(msg83)
            return ('/').join([self.out, ('').join([self.outputid, '_filtered.bam'])])

    def bam2bedgraph(self,bamfile):
	self.filelog = self.logopen()
	fnull = open(os.devnull, 'w')
	try:
	    with open( ("/").join([self.out, ('').join([self.outputid, '.bedgraph'])]), 'w') as bedgraph:
		subprocess.check_call(['bedtools','genomecov','-bg','-ibam', bamfile,'-g',self.genomesize], stderr=fnull, stdout=bedgraph)
        except traceback:
            self.filelog.write(msg82)
            sys.exit(1)
        else:
            self.filelog.write(msg83)
            return ("/").join([self.out, ('').join([self.outputid, '.bedgraph'])])

    def bedgraph2bw(self,bedgraph):
        self.filelog = self.logopen()
        fnull = open(os.devnull, 'w')
        try:
           subprocess.check_call(['bedGraphToBigWig', bedgraph, self.genomesize,self.out + ".bw"])
        except traceback:
            self.filelog.write(msg82)
            sys.exit(1)
        else:
            self.filelog.write(msg83)
            #return ("/").join([self.out, ('').join([self.outputid, '.bw'])])
            return self.out + ".bw"



    def kallistoindexvalues(self, fastaconcatenate):
        try:
            a = []
            for seq_record in SeqIO.parse(fastaconcatenate, 'fastq'):
                a.append(len(seq_record))

        except traceback:
            self.filelog.write(msg131)
            sys.exit(1)
        else:
            self.filelog.write(msg132)
            return (
             np.mean(a), np.std(a))

    def indexkallisto(self):
        self.filelog = self.logopen()
        fnull = open(os.devnull, 'w')
        try:
            if os.path.exists(self.dbname):
                self.filelog.write(msg123)
            else:
                subprocess.check_call(['kallisto', 'index', '-i', self.dbname, self.fastasequence], stderr=fnull, stdout=fnull)
        except subprocess.CalledProcessError:
            self.filelog.write(msg43)
            sys.exit(1)
        except traceback:
            self.filelog.write(msg43)
            sys.exit(1)
        else:
            self.filelog.write(msg44)
            return self.dbname

    def mappingkallistotran(self, indexkallisto,avefragmlen, stdfragmlen, reads):
        self.filelog = self.logopen()
        fnull = open(os.devnull, 'w')
        try:
	   subprocess.check_call(['kallisto', 'quant', '-i', indexkallisto,'--single','--pseudobam', '-l', str(avefragmlen), '-s', str(stdfragmlen),"-o","".join([self.outputfolder, ('/').join(["".join([self.outputid,"tr"]), ''])]),reads], stdout=fnull,stderr=fnull)
        except subprocess.CalledProcessError:
            self.filelog.write(msg60)
            sys.exit(1)
        else:
            self.filelog.write(msg61)
            return "".join([self.outputfolder, ('').join(["".join([self.outputid,"tr"]), '/pseudoalignments.bam'])])


    def mappingkallistopaired(self, indexkallisto,gtf,genomesize,read1, read2):
        self.filelog = self.logopen()
        fnull = open(os.devnull, 'w')
        try:
           subprocess.check_call(['kallisto', 'quant', '-i', indexkallisto, '-o', ('/').join([self.outputfolder, self.outputid]),'--genomebam', '--gtf', gtf, '--chromosomes', genomesize, read1, read2], stdout=fnull,stderr=fnull)
        except subprocess.CalledProcessError:
            self.filelog.write(msg60)
            sys.exit(1)
        else:
            self.filelog.write(msg61)
            return ('/').join([self.outputfolder, ('/').join([self.outputid, 'pseudoalignments.bam'])])

    def mappingkallisto(self, indexkallisto, gtf,genomesize,avefragmlen, stdfragmlen, reads):
        self.filelog = self.logopen()
        fnull = open(os.devnull, 'w')
        try:
            subprocess.check_call(['kallisto', 'quant', '-i', indexkallisto, '-o', ('/').join([self.outputfolder, self.outputid]),'--genomebam', '--gtf', gtf, '--chromosomes', genomesize, '--single', '-l', str(avefragmlen), '-s', str(stdfragmlen), reads], stdout=fnull,stderr=fnull)
        except subprocess.CalledProcessError:
            self.filelog.write(msg60)
            sys.exit(1)
        else:
            self.filelog.write(msg61)
            return ('/').join([self.outputfolder, ('/').join([self.outputid, 'pseudoalignments.bam'])])

    def mappingkallistopairedttran(self, indexkallisto, read1, read2):
        self.filelog = self.logopen()
        fnull = open(os.devnull, 'w')
        try:
           subprocess.check_call(['kallisto', 'quant', '-i', indexkallisto,'--pseudobam', '-o', ('/').join([self.outputfolder, "".join([self.outputid,"tr"])]),read1, read2], stdout=fnull,stderr=fnull)
        except subprocess.CalledProcessError:
            self.filelog.write(msg60)
            sys.exit(1)
        else:
            self.filelog.write(msg61)
            return ('/').join([self.outputfolder, ('/').join(["".join([self.outputid,"tr"]), "pseudoalignments.bam"])])

    def convert(self, samfile):
        self.filelog = self.logopen()
        try:
            pysam.view("-hbS", samfile, "-F", "4", "-@", self.thread, "-o", self.out + ".bam", catch_stdout=False)
        except traceback:
            self.filelog.write(msg82)
            sys.exit(1)
        else:
            self.filelog.write(msg83)
            return self.out + ".bam"

    def sortbam(self, bamfile):
        """
        Sort and indexing of BAM file
        :param bamfile:
        :return:
        """
        self.filelog = self.logopen()
        try:
            pysam.sort("-@", self.thread, bamfile, "-o", self.out + "_mapping.bam")
            pysam.index(self.out + "_mapping.bam")
        except subprocess.CalledProcessError:
            self.filelog.write(msg84)
            sys.exit(1)
        else:
            self.filelog.write(msg85)
            return self.out + "_mapping.bam"


    def cleantempfile(self):
        """
        Remove temporany files.
        :return:
        """
        templistfilesingle = [
         '.sam', '_filtered.sam', '.bam', '_forward.fastq', '_forward5trimmed.fastq',
         '_forward3trimmed.fastq', '_reverse3trimmed.fasta', '_forward3trimmed.fasta',
         '_reverse5trimmed.fasta', '_forward5trimmed.fasta', '_con.fasta']
        templistfilepaired = ['_forward.fastq', '_reverse.fastq', '_forward5trimmed.fastq', '_reverse5trimmed.fastq',
         '_forward3trimmed.fastq', '_reverse3trimmed.fastq', '_con.fastq', '.sam', '_filtered.sam',
         '.bam']
        if os.path.isfile(self.outputfolder + 'abundance.tsv'):
            os.remove(self.outputfolder + self.outputid + 'abundance.tsv')
        if os.path.isfile(self.outputfolder +  self.outputid + 'run_info.json'):
            os.remove(self.outputfolder +  self.outputid + 'run_info.json')
        if os.path.isfile(self.outputfolder +  self.outputid + 'abundance.h5'):
            os.remove(self.outputfolder +  self.outputid + 'abundance.h5')
	subprocess.check_call(['rm','-rf',self.outputfolder + self.outputid + "tr"])
	subprocess.check_call(['rm','-rf',self.outputfolder + self.outputid])
        if self.sequencingtype in 'Single-End':
            if self.readforwardtype in 'fastq':
                for item in templistfilesingle:
                    if os.path.isfile(self.out + item):
                        os.remove(self.out + item)

            elif self.readforwardtype in 'fasta':
                for item in templistfilesingle:
                    if os.path.isfile(self.out + item):
                        os.remove(self.out + item)

        else:
            if self.sequencingtype in 'Paired-End':
                if self.readforwardtype in 'fastq':
                    for item in templistfilepaired:
                        if os.path.isfile(self.out + item):
                            os.remove(self.out + item)

                elif self.readforwardtype in 'fasta':
                    for item in templistfilepaired:
                        if os.path.isfile(self.out + item):
                            os.remove(self.out + item)
