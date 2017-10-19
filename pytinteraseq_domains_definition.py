from pyinteraseq_inputcheck import InputCheck
from output_message import *
import pandas as pd
import sys
import subprocess
import re
import pybedtools

class DomainsDefinition(InputCheck):

    def __init__(self, optparseinstance):
        InputCheck.__init__(self, optparseinstance)
        self.df = None
        self.dfOp = None
        self.dfMM = None
        self.dflen = None
        self.dfstart = None
        self.dfMerge = None

    def hashclean(self, blastnout):
        with open(blastnout) as oldfile, open(self.out + 'blastn_nohash.tab', 'w') as newfile:
            for line in oldfile:
                if not line.startswith('#'):
                    newfile.write(line)
        return self.out + 'blastn_nohash.tab'

    def blastnfiltering(self, blastnout):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        self.df = pd.read_csv(blastnout, sep='\t', header=None,
                              names=['seq', 'chr', 'percmatch', 'length', 'mismatch', 'op', 'cstart', 'cend',
                                     'start', 'end', 'evalue', 'bitscore', 'nseq'])
        self.filelog.write(msg62 + str(len(self.df)))
        # filter open gap
        self.dfOp = self.df[(self.df.op <= 1)]
        self.filelog.write(msg63 + str(len(self.dfOp)))
        # trasform mismatch in percentage of mismatch using clone length
        self.dfOp['mismatch'] = (self.dfOp['mismatch'] / self.dfOp['length']) * 100
        # self.dfOp.loc[self.dfOp.mismatch > 0,'mismatch'] = self.dfOp['mismatch']*100 /self.dfOp['length'] #TODO
        # trasform into numeric field
        self.dfOp[['mismatch']].apply(pd.to_numeric)
        # filter on percentage of mismatch
        self.dfMM = self.dfOp[(self.dfOp['mismatch'] < 5.1)]
        self.filelog.write(msg64 + str(len(self.dfMM)))
        # lenght filtering
        self.dflen = self.dfMM[(self.dfMM.length >= 149)]
        self.filelog.write(msg65 + str(len(self.dflen)))
        # filter in start clone
        self.dfstart = self.dflen[(self.dflen.cstart <= 1)]
        # drop duplicate
        self.dfstart = self.dfstart.drop_duplicates(subset='seq', keep=False)
        self.filelog.write(msg66 + str(len(self.dfstart)))
        if self.sequencingtype == 'Paired-End':
            # split field seq in two columns
            self.dfstart['read'], self.dfstart['seqid'] = self.dfstart['seq'].str.split('_', 2).str[0:2].str
            # split into two df read1 and read2
            self.df1 = self.dfstart[(self.dfstart['read'] == 'seq1')]
            self.df2 = self.dfstart[(self.dfstart['read'] == 'seq2')]
            # merge df
            self.dfMerge = pd.merge(self.df1, self.df2, on='seqid')
            # write output
            self.dfMerge[['seq_x', 'nseq_x']].to_csv(self.out + '_p1.tab', header=None, sep='\t',
                                                     index=False)
            self.dfMerge[['seq_y', 'nseq_y']].to_csv(self.out + '_p2.tab', header=None, sep='\t',
                                                     index=False)
            self.dfMerge.to_csv(self.out + '_filtered_pairs.tab', header=None, sep='\t', index=False)
            return self.out + '_filtered_pairs.tab'
        else:
            self.dfstart.to_csv(self.out + '_filtered_single_complete.tab', header=None, sep='\t', index=False)
            self.dfstart[['seq', 'nseq']].to_csv(self.out + '_filtered_single.tab', header=None, sep='\t', index=False)
            return self.out + '_filtered_single.tab'

    def clustering(self, blastnout):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            subprocess.check_call(['pick_otus.py', '-i', blastnout, '-o', self.out + '_picked', '-s', '0.97', '-g', '10'])
        except subprocess.CalledProcessError:
            self.filelog.write(msg71)
            sys.exit(0)
        return self.out + '_picked/' + self.outputid + '_newid_otus.txt'

    def pickrepseq(self, pickotus, fasta):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            subprocess.check_call(
                ['pick_rep_set.py', '-i', pickotus, '-f',
                 fasta,
                 '-m', 'most_abundant', '-o', self.out + '_otus_most_abundant.fa'])
        except subprocess.CalledProcessError:
            self.filelog.write(msg72)
            sys.exit(0)
        else:
            self.filelog.write(msg73)
            return self.out + '_otus_most_abundant.fa'

    def pysed(self, pickrepseqoutput, id, old, new):
        with open(pickrepseqoutput) as f:
            with open(self.out+id, 'w') as g:
                for line in f:
                    g.write(re.sub(old, new, line))
            g.close()
        return self.out+id

    def bedparsing(self, blastnclonesinput):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        self.df1 = pd.read_csv(blastnclonesinput, sep="\t", header=None, names=['chr', 'start', 'end', 'clonename', 'score', 'strand'])
        self.df2 = self.df1.replace({'minus': '-', 'plus': '+'}, regex=True).sort_values('start').reset_index(drop=True)
        self.dfplus = self.df2.loc[self.df2['strand'] == '+'].reset_index(drop=True)
        self.dfminus = self.df2.loc[self.df2['strand'] == '-'].reset_index(drop=True)
        self.dfminus = self.dfminus[['chr', 'end', 'start', 'clonename', 'score', 'strand']]
        self.dfminus = self.dfminus.rename(
            columns={'chr': 'chr', 'end': 'start', 'start': 'end', 'clonesname': 'clonesname', 'score': 'score',
                     'strand': 'strand'})
        self.df1 = pd.concat([self.dfplus, self.dfminus], ignore_index=True).sort_values('start').reset_index(drop=True)
        self.df1.to_csv(self.out+'_blastclonesparsed.bed', sep="\t", header=None, index=False)
        self.filelog.write(msg74)
        return self.out+'_blastclonesparsed.bed'

    def clonescount(self, pickotusout):
        with open(self.out+'_cluster_count.txt', 'w') as f:
            subprocess.check_call(['awk', 'BEGIN{FS="\t";OFS="\t"}{print $1,NF}', pickotusout], stdout=f)
        return self.out+'_cluster_count.txt'

    def mergingcount(self, bedparsed, clonescounted):
        self.df1 = pd.read_csv(bedparsed, sep="\t", header=None, names=['chr', 'start', 'end', 'clonename', 'score', 'strand'])
        self.df2 = pd.read_csv(clonescounted, sep="\t", header=None, names=['clonename', 'count'])
        self.df3 = pd.merge(self.df1, self.df2, on='clonename')
        self.df3 = self.df3[['chr', 'start', 'end', 'clonename', 'count', 'strand']]
        self.df3.to_csv(self.out+'_blastnclonescounted.bed', sep="\t", header=None, index=False)
        return self.out+'_blastnclonescounted.bed'

    def filteringclonescount(self, mergingcountoutput, frequency):
        self.df1 = pd.read_csv(mergingcountoutput, sep="\t", header=None, names=['chr', 'start', 'end', 'clonename', 'count', 'strand'])
        self.df2 = self.df1.loc[self.df1['count'] >= int(frequency)].sort_values('start').reset_index(drop=True)
        self.df2.to_csv(self.out+'_blastnclonescountedfiltered.bed', sep="\t", header=None, index=False)
        return self.out+'_blastnclonescountedfiltered.bed'

    def pybedtoolsmerge(self, filteringclonescountoutput):
        self.input = pybedtools.example_bedtool(filteringclonescountoutput)
        self.output = self.input.merge().moveto(self.out+'_blastnclonesmerge.bed')
        return self.out+'_blastnclonesmerge.bed'

    def pybedtoolstofasta(self, pybedtoolsmergeoutput, fastqsequence):
        self.bed = pybedtools.BedTool(pybedtoolsmergeoutput)
        self.fasta = pybedtools.BedTool(fastqsequence)
        self.end = self.bed.sequence(fi=self.fasta).save_seqs(self.out + '_blastnclonesmerge.fasta')
        return self.out + '_blastnclonesmerge.fasta'

    def annotationclones(self):



