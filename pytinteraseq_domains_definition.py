from pyinteraseq_inputcheck import InputCheckDomainDefinition
from output_message import *
import pandas as pd
import sys
import subprocess
import re
import pybedtools
import traceback
import os


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
        self.fasta = None
        self.clones = None
        self.intersection = None
        self.dfplus = None
        self.dfminus = None
        self.seqix = 1
        self.id = '1'
        self.id = [(int(x) - 1) for x in self.id]
        self.header = None
        self.path_multiblastn = os.path.dirname(os.path.realpath(__file__)) + '/pyinteraseq_multblastn.py'

    def tab2fasta(self, tabular, prefixoutput):
        """
        Convert tabular fasta file (2 columns) first ID second nucleotide sequence in fasta file format
        :param tabular: File fasta in tabular format
        :param prefixoutput: prefix to append
        :return: path + idanalysis + prefix + .fasta
        """
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
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
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg100)
            sys.exit(1)
        else:
            self.filelog.write(msg101)
            return self.out + prefixoutput + '.fasta'

    def fasta2tabular(self, imp, prefix):
        """
        Function for conversion Fasta file in tabular format
        :param imp: input file in fasta format
        :param prefix: prefix to add at converted file
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            SeqIO.convert(imp, 'fasta', self.out + prefix + '.tab', 'tab')
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg96)
            sys.exit(1)
        else:
            self.filelog.write(msg97)
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
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            subprocess.check_call(['python', self.path_multiblastn,
                                   '--referencefasta', fasta,
                                   '--multifastasequence', multifasta,
                                   '--dbname', self.chromosomename,
                                   '--outputfolder', self.outputfolder,
                                   '--outputid', self.outputid + suffix,
                                   '--thread', self.thread,
                                   '--outformat', outputformat,
                                   '--log', self.outputfolder + self.outputid + ".log"],
                                  stderr=self.filelog, stdout=self.filelog)
        except subprocess.CalledProcessError:
            self.filelog.write(msg60)
            sys.exit(1)
        else:
            self.filelog.write(msg61)
            return self.out + suffix

    def clustering(self, blastnout, prefixoutput):
        """
        Perform read clustering algoritm
        :param blastnout: Table generated with filtering function
        :param prefixoutput: Prefix that will add to output
        :return: cluster file
        """
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            subprocess.check_call(['pick_otus.py', '-i', blastnout, '-o',
                                   self.out + '_picked', '-s', '0.97'])
        except subprocess.CalledProcessError:
            self.filelog.write(msg71)
            sys.exit(1)
        return self.out + '_picked/' + self.outputid + prefixoutput + '_otus.txt'

    def pickrepseq(self, pickotus, fasta):
        """
        Perform that selection of most representative sequence selection the most abundant in the cluster
        :param pickotus: Output generated by clustering function
        :param fasta: File generated with tab2fasta containing fasta sequence filtered genarated by blastn align
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            subprocess.check_call(
                ['pick_rep_set.py', '-i', pickotus, '-f',
                 fasta,
                 '-m', 'most_abundant', '-o', self.out + '_otus_most_abundant.fa'])
        except subprocess.CalledProcessError:
            self.filelog.write(msg72)
            sys.exit(1)
        else:
            self.filelog.write(msg73)
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
            self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
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
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg104)
            sys.exit(1)
        else:
            self.filelog.write(msg105)
            return self.out+'_blastclonesparsed.bed'

    def clonescount(self, pickotusout):
        """
        Function that count the dimension of cluster
        :param pickotusout:
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            with open(self.out + '_cluster_count.txt', 'w') as f:
                subprocess.check_call(['awk', 'BEGIN{FS="\t";OFS="\t"}{print $1,NF}', pickotusout], stdout=f)
        except subprocess.CalledProcessError:
            self.filelog.write(msg106)
            sys.exit(1)
        else:
            self.filelog.write(msg107)
            return self.out+'_cluster_count.txt'

    def mergingcount(self, bedparsed, clonescounted):
        """
        Function for merging Bed with analysis information to cluster clones count
        :param bedparsed:
        :param clonescounted:
        :return: BED6 files with chr, start, end, clonename, count, strand
        """
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            self.df1 = pd.read_csv(bedparsed, sep="\t", header=None,
                                   names=['chr', 'start', 'end', 'clonename', 'score', 'strand'])
            self.df2 = pd.read_csv(clonescounted, sep="\t", header=None, names=['clonename', 'count'])
            self.df3 = pd.merge(self.df1, self.df2, on='clonename')
            self.df3 = self.df3[['chr', 'start', 'end', 'clonename', 'count', 'strand']]
            self.df3.to_csv(self.out + '_blastnclonescounted.bed', sep="\t", header=None, index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg108)
            sys.exit(1)
        else:
            self.filelog.write(msg109)
            return self.out + '_blastnclonescounted.bed'

    def filteringclonescount(self, mergingcountoutput, frequency):
        """
        Function for filtering cluster with dimension less than frequency
        :param mergingcountoutput:
        :param frequency:
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            self.df1 = pd.read_csv(mergingcountoutput, sep="\t", header=None,
                                   names=['chr', 'start', 'end', 'clonename', 'count', 'strand'])
            self.df2 = self.df1.loc[self.df1['count'] >= int(frequency)].sort_values('start').reset_index(drop=True)
            self.df2.to_csv(self.out + '_blastnclonescountedfiltered.bed', sep="\t", header=None, index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg108)
            sys.exit(1)
        else:
            self.filelog.write(msg109)
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
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            self.bed = pybedtools.BedTool(pybedtoolsmergeoutput)
            self.fasta = pybedtools.BedTool(fastqsequence)
            self.end = self.bed.sequence(fi=self.fasta).save_seqs(self.out + '_blastnclonesmerge.fasta')
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg88)
            sys.exit(1)
        else:
            self.filelog.write(msg91)
            return self.out + '_blastnclonesmerge.fasta'

    def bedtoolsannotate(self, clonesformatbed, annotation):
        """
        Perform Bedtools Annotate
        :param clonesformatbed:
        :param annotation:
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            with open(self.out + '_clonesannotated.bed', 'w') as f:
                subprocess.check_call(['bedtools', 'annotate', '-i', clonesformatbed, '-files', annotation], stdout=f)
            f.close()
        except subprocess.CalledProcessError:
            self.filelog.write(msg75)
            sys.exit(0)
        else:
            self.filelog.write(msg76)
            return self.out + '_clonesannotated.bed'

    def bedtoolsannotatefiltering(self, bedtoolsannotateout, percthr):
        """

        :param bedtoolsannotateout:
        :param percthr:
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
            self.df1 = pd.read_csv(bedtoolsannotateout, sep="\t", header=None)
            self.df2 = self.df1.loc[self.df1[6] >= float(percthr)].sort_values(1).reset_index(drop=True)
            self.df2[[0, 1, 2, 3, 4, 5]].to_csv(self.out + '_clonesannotatedfiltered.bed', sep="\t",
                                                header=None, index=False)
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg88)
            sys.exit(1)
        else:
            self.filelog.write(msg77)
            return self.out + '_clonesannotatedfiltered.bed'

    def adddescription(self, clonesmerged, annotation, percthr):
        """
        Function that attach description to bed6 table
        :param clonesmerged: bed6 with clone information
        :param annotation: annotation file
        :param percthr: Overlap intersection
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
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
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg110)
            sys.exit(1)
        else:
            self.filelog.write(msg111)
            return self.out + '_clonesdescription.bed'

    def addsequence(self, outputfromdescription, outputfasta2tab):
        """

        :param outputfromdescription:
        :param outputfasta2tab:
        :return:
        """
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
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
            self.filelog.write(traceback.format_exc())
            self.filelog.write(msg112)
            sys.exit(1)
        else:
            self.filelog.write(msg113)
            return self.out + '_domaindetection_step1.tab'

    # def cleantemporaryfilesinglend(self, filedict):
    #     for key, value in filedict.iteritems():
    #         if key == 'Trimmed5single':
    #             os.remove(value)
    #         elif key == 'TabularRenamedForward':
    #             os.remove(value)
    #         elif key == 'FastaReadsForward':
    #             os.remove(value)
    #         elif key == 'FastaRenamedForward':
    #             os.remove(value)
    #         elif key == 'blastoutputnohash':
    #             os.remove(value)
    #         elif key == 'pickedreads':
    #             os.remove(value)
    #         elif key == 'blastoutputnohashfiltered':
    #             os.remove(value)
    #         elif key == 'pickedreadscleand':
    #             os.remove(value)
    #         elif key == 'blastedclones':
    #             os.remove(value)
    #         elif key == 'bedparsed':
    #             os.remove(value)
    #         elif key == 'clonesannotated':
    #             os.remove(value)
    #         elif key == 'clustercount':
    #             os.remove(value)
    #         elif key == 'clonescounted':
    #             os.remove(value)
    #         elif key == 'clonenseqfasta':
    #             os.remove(value)
    #         elif key == 'clonesmergedfasta':
    #             os.remove(value)
    #         elif key == 'clonescountedmerged':
    #             os.remove(value)
    #         elif key == 'tabwithdescription':
    #             os.remove(value)
    #         elif key == 'clonescountedfiltered':
    #             os.remove(value)
    #         elif key == 'blastoutput':
    #             os.remove(value)
    #         elif key == 'blastoutputnohashfilteredfasta':
    #             os.remove(value)
    #         shutil.rmtree(self.out + '_picked', ignore_errors=True)
    #         filelist = [f for f in os.listdir(self.outputfolder) if f.startswith(self.chromosomename)]
    #         for f in filelist:
    #             os.remove(os.path.join(self.outputfolder, f))
    #     return 0
