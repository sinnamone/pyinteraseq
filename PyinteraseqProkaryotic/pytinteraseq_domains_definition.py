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
        self.path_pickotus = os.path.dirname(os.path.realpath(__file__)) + '/pick_otus.sh'
        #self.pythoneve = "/opt/miniconda3/envs/qiime1/bin/python"


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
        self.df1[[u'seq',u'nseq']].to_csv(self.out + '_mappingoutput.tab', sep="\t",header=None,index=False)
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
                                   '--dbname', self.outputid,
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
            #"$pickotus" -i "$blastnoutput" -o "$outputfoler" -s 0.97
            subprocess.check_call(['pick_otus.py','-i',blastnout,'-o',self.out + '_picked','-s','0.97','-m','sumaclust','--threads', self.thread])
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
                ['pick_rep_set.py', '-i', pickotus, '-f',
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
	    with open(self.out + '_inputformerge.bed', 'w') as f:
		subprocess.check_call(['cut','-f','1,2,3,5',self.out + '_blastnclonescountedfiltered.bed'],stdout=f)
	    f.close()
            with open(self.out + '_inputfobigwig.bed', 'w') as f:
            	subprocess.check_call(['bedtools','merge','-c','4','-o','sum','-i',self.out + '_inputformerge.bed'],stdout=f)
            f.close()
        except traceback:
            self.filelogerrorwrite(msg108)
        else:
            self.filelogstdoutwrite(msg109)
            return self.out+'_blastnclonescountedfiltered.bed'
    
    def bigwigcreation(self):
    	"""
    	Bigwig creation
    	"""
	try:
		with open(self.outputfolder + self.namefilefasta.split('.')[0] + '.genome', "wb") as self.genome:
	    		for rec in SeqIO.parse(self.fastasequence,"fasta"):
     				self.genome.write(rec.id+"\t"+str(len(rec.seq)))
			self.genome.close()
			subprocess.check_call(['bedGraphToBigWig',self.out + '_inputfobigwig.bed',self.genome.name,self.out + '.bw'])
        except subprocess.CalledProcessError:
            self.filelogerrorwrite(msg75)
        else:
            self.filelogstdoutwrite(msg76)
            return self.out + '.bw'
			
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
                                   names=['chr', 'clonestart', 'cloneend', 'score',
                                          'start', 'end', 'geneid',
                                          'strand', 'description', 'clonelength'])
            self.df2 = pd.read_csv(outputfasta2tab, sep="\t", header=None, names=['id_tab', 'nseq'])
            self.df1['temp_id_tab'] = self.df1[['clonestart', 'cloneend']].astype(str).apply(lambda x: '-'.join(x),
                                                                                             axis=1)
            self.df1['id_tab'] = self.df1[['chr', 'temp_id_tab']].astype(str).apply(lambda x: ':'.join(x), axis=1)
            self.df1.drop('temp_id_tab', axis=1, inplace=True)
            self.df3 = pd.merge(self.df1, self.df2, on='id_tab')
            self.df3 = self.df3.loc[self.df3['geneid'] != "."].reset_index(drop=True)
            self.df4 = self.df3[['chr', 'clonestart', 'cloneend',
                                 'clonelength', 'start', 'end',
                                 'geneid', 'strand',
                                 'description', 'nseq']]
            self.df4.columns = ["Chr", "CloneStart", "CloneEnd", "CloneLength", "Start", "End", "GeneID", "Strand",
                                "Description", "NuclSeq"]
            self.df4.to_csv(self.out + '_definition.txt', sep="\t", header=True, index=False)

        except traceback:
            self.filelogerrorwrite(msg112)
        else:
            self.filelogstdoutwrite(msg113)
            return self.out + '_definition.txt'

    def cleantempfile(self):
        """
        Remove temporany files.
        :return:
        """
        db1 = str(self.out + ".nsq")
        db2 = str(self.out + ".nin")
        db3 = str(self.out + ".nhr")
        os.remove(db1)
        os.remove(db2)
        os.remove(db3)
        templistfile = ["_def_blastnfiltered.fasta", "_otus_most_abundant.fa", "_def_clean.fasta", "_def_cluster_count.txt",
                              "_def_clonestabular.tab", "_clonesdescription.bed", "_clonesannotatedfiltered.bed", "_clonesannotated.bed",
                              "_blastnclonesmerge.fasta", "_blastnclonesmerge.bed", "_blastnclonescountedfiltered.bed", "_blastnclonescounted.bed",
                              "_def_blastnclones.tab", "_def_blastclonesparsed.bed", "_blastnfiltered.fasta", "_clean.fasta",
                              "_cluster_count.txt", "_clonestabular.tab", "_blastnclones.tab", "_blastclonesparsed.bed",
                              "_mappingoutput.tab", "_inputfobigwig.bed", "_inputformerge.bed", ".genome"]
        for item in templistfile:
            if os.path.isfile(self.out + item):
                os.remove(self.out + item)
        pickedlist = ["_blastnfiltered_otus.log","_blastnfiltered_otus.txt","_blastnfiltered_clusters.uc"]
        for item in pickedlist:
            if os.path.isfile(self.out + '_picked/' + self.outputid + item):
                os.remove(self.out + '_picked/' + self.outputid + item)
        os.rmdir(self.out + '_picked/')

