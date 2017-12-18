from pyinteraseq_inputcheck import InputCheckDomainDefinition
from output_message import *
import pandas as pd
import sys
import subprocess
import re
import pybedtools
import traceback
import os
from Bio import SeqIO
import pysam


class DomainsDefinition(InputCheckDomainDefinition):

    def __init__(self, optparseinstance):
        InputCheckDomainDefinition.__init__(self, optparseinstance)
        self.df = None
        self.allPossibilities = []
        self.df1 = None
        self.df2 = None
        self.df3 = None
        self.df4 = None
        self.df5 = None
        self.input = None
        self.output = None
        self.end = None
        self.bed = None
        self.clones = None
        self.intersection = None
        self.dfplus = None
        self.dfminus = None
        self.seqix = 1
        self.id = '1'
        self.id = [(int(x) - 1) for x in self.id]
        self.header = None
        self.path_multiblastn = os.path.dirname(os.path.realpath(__file__)) + '/pyinteraseq_multblastn.py'
        self.dbname = self.outputfolder + os.path.basename(self.fastasequence.split('/')[-1]).split('.')[0]
        self.mappingoutoput = self.inputistance.mappingoutput
        self.pythoneve = "/opt/miniconda3/envs/qiime1/bin/python"
        self.outfasta = self.out + "_converted.fasta"

    def filelogstdoutwrite(self, msg):
        """
        Write information about script esecution
        :param msg:
        :return:
        """
        self.filelog = self.logopen()
        self.filelog.write(msg)

    def filelogerrorwrite(self, msg):
        """
        Write error message
        :param msg:
        :return:
        """
        self.filelog = self.logopen()
        self.filelog.write(traceback.format_exc())
        self.filelog.write(msg)
        sys.exit(1)

    # def bam2rec(self, bamsorted):
    #     """
    #     Generator to convert BAM files into Biopython SeqRecords.
    #     """
    #     bam_file = pysam.Samfile(bamsorted, "rb")
    #     for read in bam_file:
    #         seq = Seq.Seq(read.seq)
    #         if read.is_reverse:
    #             seq = seq.reverse_complement()
    #         rec = SeqRecord.SeqRecord(seq, read.qname, "", "")
    #         yield rec
    #
    def bam2fasta(self, bamfile):



    def fastq2fasta(self, fastq, nameid):
        """
        Covert Fastq in Fasta format
        :param fastq: input Fastq file
        :param nameid: output name
        :return:
        """
        self.filelog = self.logopen()
        try:
            SeqIO.convert(fastq, 'fastq', self.out + nameid + '.fasta', 'fasta')
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg94)
            sys.exit(1)
        else:
            self.filelog.write(msg95)
            return self.out + nameid + '.fasta'

    def mappingoutput2tabular(self, tabularoutput):
        """

        :param tabularoutput:
        :return:
        """
        self.df1 = pd.read_csv(tabularoutput, sep ="\t", header=None, names=[u'seq', u'chr', u'percmatch', u'length',
                                                                             u'mismatch', u'op', u'cstart',
                                                                             u'cend', u'start', u'end', u'evalue',
                                                                             u'bitscore', u'nseq'])
        self.df1 = self.df1.drop_duplicates(subset='seq', keep=False)
        self.df1[[u'seq', u'nseq']].to_csv(self.out + '_mappingoutput.tab', sep="\t",header=None,index=False)
        return self.out + '_mappingoutput.tab'

    def tab2fasta(self, tabular, prefixoutput):
        """
        Convert tabular fasta file (2 columns) first ID second nucleotide sequence in fasta file format
        :param tabular: File fasta in tabular format
        :param prefixoutput: prefix to append
        :return: path + idanalysis + prefix + .fasta
        """
        try:
            with open(tabular, 'r') as f:
                with open(self.out + prefixoutput + '.fasta', 'w') as f_out:
                    for line in f:
                        line = line.strip().split('\t')
                        self.header = '>' + '_'.join([line[i] for i in self.id])
                        f_out.write(self.header + '\n')
                        f_out.write(line[self.seqix] + '\n')
                f_out.close()
        except traceback:
            self.filelogerrorwrite(msg100)
        else:
            self.filelogstdoutwrite(msg101)
            return self.out + prefixoutput + '.fasta'

    def fasta2tabular(self, imp, prefix):
        """
        Function for conversion Fasta file in tabular format
        :param imp: input file in fasta format
        :param prefix: prefix to add at converted file
        :return:
        """
        try:
            SeqIO.convert(imp, 'fasta', self.out + prefix + '.tab', 'tab')
        except traceback:
            self.filelogerrorwrite(msg96)
        else:
            self.filelogstdoutwrite(msg97)
            return self.out + prefix + '.tab'

    def callmultiblastn(self, fasta, multifasta, outputformat, suffix):
        """
        Function to call multithread blastn,
        :param fasta: Reference or fasta that will be use for makeblastdb
        :param multifasta: Input multifasta file
        :param outputformat: Blast output format
        :param suffix: String added to outputfile
        :return: blastn output
        """
        self.filelog = self.logopen()
        fnull = open(os.devnull, 'w')
        try:
            subprocess.check_call(['python', self.path_multiblastn,
                                   '--referencefasta', fasta,
                                   '--multifastasequence', multifasta,
                                   '--dbname', self.genename,
                                   '--outputfolder', self.outputfolder,
                                   '--outputid', self.outputid + suffix,
                                   '--thread', self.thread,
                                   '--outformat', outputformat])
        except subprocess.CalledProcessError:
            self.filelog.write(msg60)
            sys.exit(1)
        else:
            self.filelog.write(msg61)
            self.filelog.close()
            return self.out + suffix

    def clustering(self, blastnout, prefixoutput):
        """
        Perform read clustering algoritm
        :param blastnout: Table generated with filtering function
        :param prefixoutput: Prefix that will add to output
        :return: cluster file
        """
        try:
            subprocess.check_call([self.pythoneve, self.pick_otus, '-i', blastnout, '-o',
                                   self.out + '_picked', '-s', '0.97'])
        except subprocess.CalledProcessError:
            self.filelogerrorwrite(msg71)
        else:
            return self.out + '_picked/' + self.outputid + prefixoutput + '_otus.txt'

    def pickrepseq(self, pickotus, fasta):
        """
        Perform that selection of most representative sequence selection the most abundant in the cluster
        :param pickotus: Output generated by clustering function
        :param fasta: File generated with tab2fasta containing fasta sequence filtered genarated by blastn align
        :return:
        """
        try:
            subprocess.check_call(
                [self.pythoneve,self.pick_rep_set, '-i', pickotus, '-f',
                 fasta,
                 '-m', 'most_abundant', '-o', self.out + '_otus_most_abundant.fa'])
        except subprocess.CalledProcessError:
            self.filelogerrorwrite(msg72)
        else:
            self.filelogstdoutwrite(msg73)
            return self.out + '_otus_most_abundant.fa'

    def pysed(self, pickrepseqoutput, idx, old, new):
        """
        Sed function
        :param pickrepseqoutput: Output generated by pickrepseq function
        :param idx: prefix that will be add to output file name
        :param old: string that will be substituted
        :param new: new string
        :return:
        """
        with open(pickrepseqoutput) as f:
            with open(self.out+idx, 'w') as g:
                for line in f:
                    g.write(re.sub(old, new, line))
            g.close()
        return self.out + idx

    def bedparsing(self, blastnclonesinput):
        """

        :param blastnclonesinput:
        :return:
        """
        try:
            self.df1 = pd.read_csv(blastnclonesinput, sep="\t", header=None,
                                   names=['chr', 'start', 'end', 'clonename', 'score', 'strand'])
            self.df2 = self.df1.replace({'minus': '-', 'plus': '+'}, regex=True).sort_values('start').reset_index(
                drop=True)
            self.dfplus = self.df2.loc[self.df2['strand'] == '+'].reset_index(drop=True)
            self.dfminus = self.df2.loc[self.df2['strand'] == '-'].reset_index(drop=True)
            self.dfminus = self.dfminus[['chr', 'end', 'start', 'clonename', 'score', 'strand']]
            self.dfminus = self.dfminus.rename(
                columns={'chr': 'chr', 'end': 'start', 'start': 'end', 'clonesname': 'clonesname', 'score': 'score',
                         'strand': 'strand'})
            self.df1 = pd.concat([self.dfplus, self.dfminus],
                                 ignore_index=True).sort_values('start').reset_index(drop=True)
            self.df1.to_csv(self.out + '_blastclonesparsed.bed', sep="\t", header=None, index=False)
        except traceback:
            self.filelogerrorwrite(msg104)
        else:
            self.filelogstdoutwrite(msg105)
            return self.out+'_blastclonesparsed.bed'

    def clonescount(self, pickotusout):
        """
        Function that count the dimension of cluster
        :param pickotusout:
        :return:
        """
        try:
            with open(self.out + '_cluster_count.txt', 'w') as f:
                subprocess.check_call(['awk', 'BEGIN{FS="\t";OFS="\t"}{print $1,NF}', pickotusout], stdout=f)
        except subprocess.CalledProcessError:
            self.filelogerrorwrite(msg106)
        else:
            self.filelogstdoutwrite(msg107)
            return self.out+'_cluster_count.txt'

    def mergingcount(self, bedparsed, clonescounted):
        """
        Function for merging Bed with analysis information to cluster clones count
        :param bedparsed:
        :param clonescounted:
        :return: BED6 files with chr, start, end, clonename, count, strand
        """
        try:
            self.df1 = pd.read_csv(bedparsed, sep="\t", header=None,
                                   names=['chr', 'start', 'end', 'clonename', 'score', 'strand'])
            self.df2 = pd.read_csv(clonescounted, sep="\t", header=None, names=['clonename', 'count'])
            self.df3 = pd.merge(self.df1, self.df2, on='clonename')
            self.df3 = self.df3[['chr', 'start', 'end', 'clonename', 'count', 'strand']]
            self.df3.to_csv(self.out + '_blastnclonescounted.bed', sep="\t", header=None, index=False)
        except traceback:
            self.filelogerrorwrite(msg108)
        else:
            self.filelogstdoutwrite(msg109)
            return self.out + '_blastnclonescounted.bed'

    def filteringclonescount(self, mergingcountoutput, frequency):
        """
        Function for filtering cluster with dimension less than frequency
        :param mergingcountoutput:
        :param frequency:
        :return:
        """
        try:
            self.df1 = pd.read_csv(mergingcountoutput, sep="\t", header=None,
                                   names=['chr', 'start', 'end', 'clonename', 'count', 'strand'])
            self.df2 = self.df1.loc[self.df1['count'] >= int(frequency)].sort_values('start').reset_index(drop=True)
            self.df2 = self.df2.sort_values(['chr', 'start'], ascending=[True, True])
            self.df2 = self.df2.drop_duplicates(subset='start', keep="last")
            self.df2 = self.df2.drop_duplicates(subset='end', keep="last")
            self.df2.to_csv(self.out + '_blastnclonescountedfiltered.bed', sep="\t", header=None, index=False)
        except traceback:
            self.filelogerrorwrite(msg108)
        else:
            self.filelogstdoutwrite(msg109)
            return self.out+'_blastnclonescountedfiltered.bed'

    def pybedtoolsmerge(self, filteringclonescountoutput):
        self.input = pybedtools.example_bedtool(filteringclonescountoutput)
        self.output = self.input.merge().moveto(self.out+'_blastnclonesmerge.bed')
        return self.out+'_blastnclonesmerge.bed'

    def pybedtoolstofasta(self, pybedtoolsmergeoutput, fastqsequence):
        """
        Function that produce Fasta file from Bed
        :param pybedtoolsmergeoutput: Output of function pybedtoolsmerge
        :param fastqsequence: Input fasta file
        :return:
        """
        try:
            self.bed = pybedtools.BedTool(pybedtoolsmergeoutput)
            self.fasta = pybedtools.BedTool(fastqsequence)
            self.end = self.bed.sequence(fi=self.fasta).save_seqs(self.out + '_blastnclonesmerge.fasta')
        except traceback:
            self.filelogerrorwrite(msg88)
        else:
            self.filelogstdoutwrite(msg91)
            return self.out + '_blastnclonesmerge.fasta'

    def bedtoolsannotate(self, clonesformatbed, annotation):
        """
        Perform Bedtools Annotate
        :param clonesformatbed:
        :param annotation:
        :return:
        """
        try:
            with open(self.out + '_clonesannotated.bed', 'w') as f:
                subprocess.check_call(['bedtools', 'annotate', '-i', clonesformatbed, '-files', annotation], stdout=f)
            f.close()
        except subprocess.CalledProcessError:
            self.filelogerrorwrite(msg75)
        else:
            self.filelogstdoutwrite(msg76)
            return self.out + '_clonesannotated.bed'

    def bedtoolsannotatefiltering(self, bedtoolsannotateout, percthr):
        """

        :param bedtoolsannotateout:
        :param percthr:
        :return:
        """
        try:
            self.df1 = pd.read_csv(bedtoolsannotateout, sep="\t", header=None)
            self.df2 = self.df1.loc[self.df1[6] >= float(percthr)].sort_values(1).reset_index(drop=True)
            self.df2[[0, 1, 2, 3, 4, 5]].to_csv(self.out + '_clonesannotatedfiltered.bed', sep="\t",
                                                header=None, index=False)
        except traceback:
            self.filelogerrorwrite(msg88)
        else:
            self.filelogstdoutwrite(msg77)
            return self.out + '_clonesannotatedfiltered.bed'

    def adddescription(self, clonesmerged, annotation, percthr):
        """
        Function that attach description to bed6 table
        :param clonesmerged: bed6 with clone information
        :param annotation: annotation file
        :param percthr: Overlap intersection
        :return:
        """
        try:
            self.clones = pybedtools.BedTool(clonesmerged)
            self.annotation = pybedtools.BedTool(annotation)
            self.intersection = self.clones.intersect(self.annotation, wao=True, f=float(percthr))
            self.df1 = pd.read_table(self.intersection.fn,
                                     names=['chr', 'clonestart', 'cloneend',
                                            'chr2', 'start', 'end', 'geneid',
                                            'cog', 'strand', 'genename', 'description', 'clonelength'])
            self.df2 = self.df1.loc[self.df1['clonelength'] != int(0)].sort_values('clonestart').reset_index(drop=True)
            self.df2[['chr', 'clonestart', 'cloneend',
                      'clonelength', 'start', 'end', 'geneid',
                      'strand', 'genename', 'description']].to_csv(self.out + '_clonesdescription.bed', sep="\t",
                                                                   header=None, index=False)
        except traceback:
            self.filelogerrorwrite(msg110)
        else:
            self.filelogstdoutwrite(msg111)
            return self.out + '_clonesdescription.bed'

    def addsequence(self, outputfromdescription, outputfasta2tab):
        """

        :param outputfromdescription:
        :param outputfasta2tab:
        :return:
        """
        try:
            # ex clones description sequence
            self.df1 = pd.read_csv(outputfromdescription,
                                   sep="\t", header=None,
                                   names=['chr', 'clonestart', 'cloneend', 'clonelength',
                                          'start', 'end', 'geneid',
                                          'strand', 'genename', 'description'])
            self.df2 = pd.read_csv(outputfasta2tab, sep="\t", header=None, names=['id_tab', 'nseq'])
            self.df1['temp_id_tab'] = self.df1[['clonestart', 'cloneend']].astype(str).apply(lambda x: '-'.join(x),
                                                                                             axis=1)
            self.df1['id_tab'] = self.df1[['chr', 'temp_id_tab']].astype(str).apply(lambda x: ':'.join(x), axis=1)
            self.df1.drop('temp_id_tab', axis=1, inplace=True)
            self.df3 = pd.merge(self.df1, self.df2, on='id_tab')
            self.df3[['chr', 'clonestart', 'cloneend',
                      'clonelength', 'start', 'end',
                      'geneid', 'strand', 'genename',
                      'description', 'nseq']].to_csv(self.out + '_domaindetection_step1.tab', sep="\t",
                                                     header=None, index=False)
        except traceback:
            self.filelogerrorwrite(msg112)
        else:
            self.filelogstdoutwrite(msg113)
            return self.out + '_domaindetection_step1.tab'

    def cleantempfile(self):
        """
        Remove temporany files.
        :return:
        """
        db1 = str(self.dbname + ".nsq")
        db2 = str(self.dbname + ".nin")
        db3 = str(self.dbname + ".nhr")
        os.remove(db1)
        os.remove(db2)
        os.remove(db3)
        templistfilesingle = ["_def_blastnfiltered.fasta", "_otus_most_abundant.fa", "_def_clean.fasta", "_def_cluster_count.txt",
                              "_def_clonestabular.tab", "_clonesdescription.bed", "_clonesannotatedfiltered.bed", "_clonesannotated.bed",
                              "_blastnclonesmerge.fasta", "_blastnclonesmerge.bed", "_blastnclonescountedfiltered.bed", "_blastnclonescounted.bed",
                              "_def_blastnclones.tab", "_def_blastclonesparsed.bed", "_blastnfiltered.fasta", "_clean.fasta",
                              "_cluster_count.txt", "_clonestabular.tab", "_blastnclones.tab", "_blastclonesparsed.bed"]
        for item in templistfilesingle:
            if os.path.isfile(self.out + item):
                os.remove(self.out + item)





# import os
# from os import listdir
# import subprocess
# import sys
# import time
# import pandas as pd
# import numpy as np
# import optparse
#
# parser = optparse.OptionParser(usage='%prog [options] \n'
#                                      '\n -a or --bam [ Input file BAM ] '
#                                      '\n -b or --bed [ Annotation BED file ] \n'
#                                      '\n -c or --genome [ Genome format file ] '
#                                      '\n -p or --path [ Path where write output file ] '
#                                      '\n -t or --threshold [ Percentile threshold ] '
#                                      '\n -o or --output [ Output ID ] '
#                                      '\n -h [ show this help message and exit ]',
#                                description='Domain Detection script. This script computes the identification of domain/epitopes starting from BAM file. Annotation and genome file are mandatory. Simone Puccio',
#                                prog='domain_detection.py',
#                                version='1.0'
#                                )
#
# parser.add_option('-a', '--bam', action='store', dest='bam', default=False, help='Input file BAM')
# parser.add_option('-b', '--bed', action='store', dest='bed', default=False, help='Annotation BED file')
# parser.add_option('-c', '--genome', action='store', dest='genome', default=False, help='Genome format file')
# parser.add_option('-p', '--path', action='store', dest='path', default=False, help='Path where write output file')
# parser.add_option('-t', '--threshold', action='store', dest='threshold', type='int', default=30,
#                   help='Percentile threshold')
# parser.add_option('-o', '--output', action='store', dest='output', default=False, help='Output ID')
# options, args = parser.parse_args()
#
#
# class DomainsDefinition(object):
#
#     def __init__(self, optparseinstance):
#         self.inputistance = optparseinstancea
#
#     def focus_calculation(file1, file2, path_out):
#         dfA = pd.read_csv(file1, index_col=False, header=None, sep='\t')
#         dfA.columns = ["chr", "num_read"]
#         dfB = pd.read_csv(file2, index_col=False, header=None, sep='\t')
#         dfB.columns = ["chr", "max"]
#         dfC = pd.merge(dfA, dfB, on='chr')
#         newcol = np.divide(dfC["max"], dfC["num_read"], dtype='float')
#         dfC["focus"] = newcol
#         dfC.to_csv(path_out + options.output + 'temp_6.bed', header=None, sep='\t', index=False)
#
#     def percentile(file, perc, path_out):
#         # dataframe creation
#         dfA = pd.read_csv(file, index_col=False, header=None, sep='\t')
#         # labeling columns
#         dfA.columns = ['a1', 'a2', 'a3']
#         # list with values of depth
#         p = dfA.a3
#         # percentile function, value int perc is the threshold
#         v = np.percentile(p, perc)
#         dfD = dfA.loc[lambda df: df.a3 > v, :]
#         dfD.to_csv(path_out + options.output + 'temp_8.bed', header=None,
#                    sep='\t', index=False)
#
#     def merge_intervals_focus(file1, file2, file3, path_out):
#         dfA = pd.read_csv(file2, index_col=False, header=None, sep='\t')
#         dfA.columns = ['chr', 'start_int', 'end_int', 'mean_depth']
#         dfB = pd.read_csv(file1, index_col=False, header=None, sep='\t')
#         dfB.columns = ["chr", "num_read", "max", "focus"]
#         dfC = pd.read_csv(file3, index_col=False, header=None, sep='\t')
#         dfC.columns = ["chr", "start", "end", "gene_id", "point", "strand", "description"]
#         dfD = pd.merge(dfA, dfB, on='chr')
#         dfE = pd.merge(dfD, dfC, on='chr')
#         dfE.to_csv(path_out + options.output + 'domain_detection.bed', header=None,
#                    sep='\t', index=False)
#
#     def groupbyreads(file):
#         df = pd.read_csv(file, index_col=False, header=None, sep='\t')
#         df.columns = ['a', 'b', 'c', 'd', 'e', 'f']
#         df1 = df.groupby(["a"], as_index=False)["d"].count()
#         df1.to_csv(options.path + options.output + 'temp_4.bed', header=None, sep='\t', index=False)
#
#     def groupbydepth(file):
#         df = pd.read_csv(file, index_col=False, header=None, sep='\t')
#         df.columns = ['a', 'b', 'c', 'd', 'e']
#         df1 = df.groupby(["a"], as_index=False)["b"].max()
#         df1.to_csv(options.path + options.output + 'temp_5.bed', header=None, sep='\t', index=False)
#
#     def clean(dest):
#         try:
#             suffixdelete = ["temp_1.bed", "temp_2.bed", "temp_3.bed", "temp_4.bed", "temp_5.bed", "temp_6.bed",
#                             "temp_7.bed", "temp_8.bed", "temp_9.bed", "temp_10.bed"]
#             for file in os.listdir(dest):
#                 for i in range(len(suffixdelete)):
#                     if file.endswith(suffixdelete[i]):
#                         os.remove(os.path.join(dest, file))
#         except ValueError:
#             sys.stdout.write('Error. Removing temporary files. Exit\n')
#             sys.exit(0)
#         else:
#             sys.stdout.write('Removing temporary files. Complete.\n')
#
#     def groupby(file):
#         df = pd.read_csv(file, index_col=False, header=None, sep='\t')
#         df.columns = ['a', 'b', 'c']
#         df1 = df.groupby(['a'], as_index=False)['b'].min()
#         df2 = df.groupby(["a"], as_index=False)["b"].max()
#         df3 = df.groupby(["a"], as_index=False)["c"].mean()
#         dfA = pd.merge(df1, df2, on='a')
#         dfB = pd.merge(dfA, df3, on='a')
#         # print dfB.head()
#         dfB.to_csv(options.path + options.output + 'temp_7.bed', header=None, sep='\t', index=False)
#
#     def groupbytwo(file):
#         df = pd.read_csv(file, index_col=False, header=None, sep='\t')
#         df.columns = ['a', 'b', 'c']
#         df1 = df.groupby(['a'], as_index=False)['b'].min()
#         df2 = df.groupby(["a"], as_index=False)["b"].max()
#         df3 = df.groupby(["a"], as_index=False)["c"].mean()
#         dfA = pd.merge(df1, df2, on='a')
#         dfB = pd.merge(dfA, df3, on='a')
#         dfB.to_csv(options.path + options.output + 'temp_9.bed', header=None, sep='\t', index=False)
#
#
# if __name__ == '__main__' :
#     try:
#         with open(options.path + options.output + 'temp_1.bed', 'w') as a:
#             subprocess.check_call([bedtools + 'genomeCoverageBed', '-dz', '-ibam', options.bam, '-g', options.genome],
#                                   stdout=a)
#     except subprocess.CalledProcessError:
#         sys.stdout.write('Error in computing depth-of-coverage. Exit\n')
#         sys.exit(0)
#     else:
#         sys.stdout.write('Computing the depth-of-coverage complete.\n')
#
#     try:
#         with open(options.path + options.output + 'temp_2.bed', 'w') as b:
#             subprocess.check_call([bedtools + 'genomeCoverageBed', '-ibam', options.bam, '-g', options.genome],
#                                   stdout=b)
#     except subprocess.CalledProcessError:
#         sys.stdout.write('Error in computing breadth-of-coverage. Exit\n')
#         sys.exit(0)
#     else:
#         sys.stdout.write('Computing the breadth-of-coverage complete.\n')
#
#     try:
#         with open(options.path + options.output + 'temp_3.bed', 'w') as c:
#             subprocess.check_call([bedtools + 'bamToBed', '-i', options.bam], stdout=c)
#     except subprocess.CalledProcessError:
#         sys.stdout.write('Error in conversion BAM to BED. Exit\n')
#         sys.exit(0)
#     else:
#         sys.stdout.write('Computing the conversion BAM to BED complete.\n')
#
#     try:
#         groupbyreads(c.name)
#     except subprocess.CalledProcessError:
#         sys.stdout.write('Error. Grouping number of reads for each intervals failed.Exit\n')
#         sys.exit(0)
#     else:
#         sys.stdout.write('Grouping number of reads for each intervals.\n')
#
#     try:
#         groupbydepth(b.name)
#     except subprocess.CalledProcessError:
#         sys.stdout.write('Error. Grouping number of reads for each intervals failed.Exit\n')
#         sys.exit(0)
#     else:
#         sys.stdout.write('Grouping number of reads for each intervals.\n')
#
#     try:
#         focus_calculation(options.path + options.output + 'temp_4.bed', options.path + options.output + 'temp_5.bed',
#                           options.path)
#     except subprocess.CalledProcessError:
#         sys.stdout.write('Error. During focus index calculation. Exit\n')
#         sys.exit(0)
#     else:
#         sys.stdout.write('Focus index complete.\n')
#
#     try:
#         groupby(a.name)
#     except subprocess.CalledProcessError:
#         sys.stdout.write('Error. Parsing file for focus retrievement. Exit\n')
#         sys.exit(0)
#     else:
#         sys.stdout.write('Percentile calculation complete.\n')
#
#     # definizione intervalliii
#     try:
#         percentile(options.path + options.output + 'temp_1.bed', options.threshold, options.path)
#     except subprocess.CalledProcessError:
#         sys.stdout.write('Error. Percentile calculation failed. Exit\n')
#         sys.exit(0)
#     else:
#         sys.stdout.write('Percentile calculation complete.\n')
#
#     try:
#         groupbytwo(options.path + options.output + 'temp_8.bed')
#     except subprocess.CalledProcessError:
#         sys.stdout.write('Error. Parsing file for focus retrievement. Exit\n')
#         sys.exit(0)
#     else:
#         sys.stdout.write('Percentile calculation complete.\n')
#
#     try:
#         with open(options.path + options.output + 'temp_10.bed', 'w') as h:
#             subprocess.check_call(['awk',
#                                    'BEGIN{FS="\t";OFS="\t"}{split($1,a,"_");if($3-$2>"30" || $3-$2<"1200")print $1,($2+a[3]),($3+a[3]),$4}',
#                                    options.path + options.output + 'temp_9.bed'],
#                                   stdout=h)
#     except subprocess.CalledProcessError:
#         sys.stdout.write('Error parsing columns. Exit\n')
#         sys.exit(0)
#     else:
#         sys.stdout.write('Columns parsed complete\n')
#
#     try:
#         merge_intervals_focus(options.path + options.output + 'temp_6.bed',
#                               options.path + options.output + 'temp_10.bed', options.bed, options.path)
#     except subprocess.CalledProcessError:
#         sys.stdout.write('Error in merging intervals. Exit\n')
#         sys.exit(0)
#     else:
#         sys.stdout.write('Intervals merging complete\n')
#
#     try:
#         clean(options.path)
#     except Exception:
#         sys.stdout.write('Error. Delete temp file failed. Exit\n')
#         sys.exit(0)
#     else:
#         sys.stdout.write('Enrichement detection complete.\n')
#
#
