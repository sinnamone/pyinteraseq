from pyinteraseq_inputcheck import InputCheck
from pyinteraseq_trimming import TrimmingPaired
from output_message import *
import pandas as pd
import sys
import subprocess
import re
import pybedtools
import os
from Bio import SeqIO
from Bio.Alphabet import ProteinAlphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import warnings


class DomainsDefinition(InputCheck):

    def __init__(self, optparseinstance):
        InputCheck.__init__(self, optparseinstance)
        warnings.filterwarnings("ignore")
        self.df = None
        self.dfOp = None
        self.dfMM = None
        self.dflen = None
        self.dfstart = None
        self.dfMerge = None
        self.allPossibilities = []
        self.df1 = None
        self.df2 = None
        self.df3 = None
        self.dfplus = None
        self.dfminus = None
        self.input = None
        self.output = None
        self.end = None
        self.bed = None
        self.fasta = None
        self.clones = None
        self.intersection = None

    def hashclean(self, blastnout, prefix):
        with open(blastnout) as oldfile, open(self.out + prefix + '.tab', 'w') as newfile:
            for line in oldfile:
                if not line.startswith('#'):
                    newfile.write(line)
        return self.out + prefix + '.tab'

    def blastnfiltering(self, blastnout):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            self.df = pd.read_csv(blastnout, sep='\t', header=None,
                                  names=['seq', 'chr', 'percmatch', 'length', 'mismatch', 'op', 'cstart', 'cend',
                                         'start', 'end', 'evalue', 'bitscore', 'nseq'])
            self.filelog.write(msg62 + str(len(self.df)))
            # filter open gap
            self.dfOp = self.df[(self.df.op <= 1)]
            self.filelog.write(msg63 + str(len(self.dfOp)))
            # trasform mismatch in percentage of mismatch using clone length
            self.dfOp['pmismatch'] = (self.dfOp.mismatch.div(self.dfOp.length).mul(100))
            # trasform into numeric field
            self.dfOp[['pmismatch']].apply(pd.to_numeric)
            # filter on percentage of mismatch
            self.dfMM = self.dfOp[(self.dfOp['pmismatch'] < 5.1)]
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
                self.dfstart['read'], self.dfstart['seqid'] = self.dfstart['seq'].str.split(':', 2).str[0:2].str
                # split into two df read1 and read2
                self.df1 = self.dfstart[(self.dfstart['read'] == 'seq1')]
                self.df2 = self.dfstart[(self.dfstart['read'] == 'seq2')]
                self.df1['nread'] = self.df1['seq'].str.split(':', 2).str[2]
                self.df2['nread'] = self.df2['seq'].str.split(':', 2).str[2]
                # merge df
                self.dfMerge = pd.merge(self.df1, self.df2, on='nread')
                # write output
                self.dfForw = self.dfMerge[['seq_x', 'nseq_x']]
                self.dfRev = self.dfMerge[['seq_y', 'nseq_y']]
                self.dfForw = self.dfForw.rename(columns={'seq_x': 'seq', 'nseq_x': 'nseq'})
                self.dfRev = self.dfRev.rename(columns={'seq_y': 'seq', 'nseq_y': 'nseq'})
                self.dfMerge2 = self.dfForw.append(self.dfRev, ignore_index=True)
                self.dfMerge2[['seq', 'nseq']].to_csv(self.out + '_filtered_paired.tab', header=None, sep='\t',
                                                      index=False)
                self.dfMerge.to_csv(self.out + '_filtered_paired_complete.tab', header=None, sep='\t', index=False)
                return self.out + '_filtered_paired.tab'
            elif self.sequencingtype == 'Single-End':
                self.dfstart.to_csv(self.out + '_filtered_single_complete.tab', header=None, sep='\t', index=False)
                self.dfstart[['seq', 'nseq']].to_csv(self.out + '_filtered_single.tab', header=None, sep='\t', index=False)
                return self.out + '_filtered_single.tab'
        except Warning:
            self.filelog.write('\nWarning')
        else:
            self.filelog.write('\nOk')

    def clustering(self, blastnout, idx):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            subprocess.check_call(['pick_otus.py', '-i', blastnout, '-o',
                                   self.out + '_picked', '-s', '0.97', '-g', '10'])
        except subprocess.CalledProcessError:
            self.filelog.write(msg71)
            sys.exit(0)
        return self.out + '_picked/' + self.outputid + idx + '_otus.txt'

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

    def pysed(self, pickrepseqoutput, idx, old, new):
        with open(pickrepseqoutput) as f:
            with open(self.out+idx, 'w') as g:
                for line in f:
                    g.write(re.sub(old, new, line))
            g.close()
        return self.out + idx

    def bedparsing(self, blastnclonesinput):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        self.df1 = pd.read_csv(blastnclonesinput, sep="\t", header=None,
                               names=['chr', 'start', 'end', 'clonename', 'score', 'strand'])
        self.df2 = self.df1.replace({'minus': '-', 'plus': '+'}, regex=True).sort_values('start').reset_index(drop=True)
        self.dfplus = self.df2.loc[self.df2['strand'] == '+'].reset_index(drop=True)
        self.dfminus = self.df2.loc[self.df2['strand'] == '-'].reset_index(drop=True)
        self.dfminus = self.dfminus[['chr', 'end', 'start', 'clonename', 'score', 'strand']]
        self.dfminus = self.dfminus.rename(
            columns={'chr': 'chr', 'end': 'start', 'start': 'end', 'clonesname': 'clonesname', 'score': 'score',
                     'strand': 'strand'})
        self.df1 = pd.concat([self.dfplus, self.dfminus],
                             ignore_index=True).sort_values('start').reset_index(drop=True)
        self.df1.to_csv(self.out+'_blastclonesparsed.bed', sep="\t", header=None, index=False)
        self.filelog.write(msg74)
        return self.out+'_blastclonesparsed.bed'

    def clonescount(self, pickotusout):
        with open(self.out+'_cluster_count.txt', 'w') as f:
            subprocess.check_call(['awk', 'BEGIN{FS="\t";OFS="\t"}{print $1,NF}', pickotusout], stdout=f)
        return self.out+'_cluster_count.txt'

    def mergingcount(self, bedparsed, clonescounted):
        self.df1 = pd.read_csv(bedparsed, sep="\t", header=None,
                               names=['chr', 'start', 'end', 'clonename', 'score', 'strand'])
        self.df2 = pd.read_csv(clonescounted, sep="\t", header=None, names=['clonename', 'count'])
        self.df3 = pd.merge(self.df1, self.df2, on='clonename')
        self.df3 = self.df3[['chr', 'start', 'end', 'clonename', 'count', 'strand']]
        self.df3.to_csv(self.out+'_blastnclonescounted.bed', sep="\t", header=None, index=False)
        return self.out+'_blastnclonescounted.bed'

    def filteringclonescount(self, mergingcountoutput, frequency):
        self.df1 = pd.read_csv(mergingcountoutput, sep="\t", header=None,
                               names=['chr', 'start', 'end', 'clonename', 'count', 'strand'])
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

    def bedtoolsannotate(self, clonesformatbed, annotation):
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
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        self.df1 = pd.read_csv(bedtoolsannotateout, sep="\t", header=None)
        self.df2 = self.df1.loc[self.df1[6] >= float(percthr)].sort_values(1).reset_index(drop=True)
        self.df2[[0, 1, 2, 3, 4, 5]].to_csv(self.out + '_clonesannotatedfiltered.bed', sep="\t",
                                            header=None, index=False)
        self.filelog.write(msg77)
        return self.out + '_clonesannotatedfiltered.bed'

    def adddescription(self, clonesmerged, annotation, percthr):
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
        except IOError:
            self.filelog.write(msg78)
            sys.exit(0)
        else:
            self.filelog.write(msg79)
            return self.out + '_clonesdescription.bed'

    def addsequence(self, outputfromdescription, outputfasta2tab):
        # ex clones description sequence
        self.df1 = pd.read_csv(outputfromdescription,
                               sep="\t", header=None,
                               names=['chr', 'clonestart', 'cloneend', 'clonelength',
                                      'start', 'end', 'geneid',
                                      'strand', 'genename', 'description'])
        self.df2 = pd.read_csv(outputfasta2tab, sep="\t", header=None, names=['id_tab', 'nseq'])
        self.df1['temp_id_tab'] = self.df1[['clonestart', 'cloneend']].astype(str).apply(lambda x: '-'.join(x), axis=1)
        self.df1['id_tab'] = self.df1[['chr', 'temp_id_tab']].astype(str).apply(lambda x: ':'.join(x), axis=1)
        self.df1.drop('temp_id_tab', axis=1, inplace=True)
        self.df3 = pd.merge(self.df1, self.df2, on='id_tab')
        self.df3[['chr', 'clonestart', 'cloneend',
                  'clonelength', 'start', 'end',
                  'geneid', 'strand', 'genename',
                  'description', 'nseq']].to_csv(self.out + '_domaindetection_step1.tab', sep="\t",
                                                 header=None, index=False)
        return self.out + '_domaindetection_step1.tab'

    def translatednaframes(self, seq, inputfile):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        input_handle = open(inputfile, "r")
        output_handle = open(self.out + '_allframes.tab', "a+")
        allpossibilities = []
        for frame in range(3):
            trans = str(seq.seq[frame:].translate(11))
            allpossibilities.append(trans)
        for frame in range(3):
            trans = str(seq.seq.reverse_complement()[frame:].translate(11))
            allpossibilities.append(trans)
        i = 0
        for currentFrame in allpossibilities:
            i = i + 1
            currentprotein = Seq(currentFrame, alphabet=ProteinAlphabet)

            currentproteinrecord = SeqRecord(currentprotein, seq.name)
            currentproteinrecord.id = currentproteinrecord.id + "." + str(i)
            currentproteinrecord.description = seq.description + "; frame " + str(i)
            SeqIO.write(currentproteinrecord, output_handle, "tab")
        input_handle.close()
        output_handle.close()
        self.filelog.write(msg83)
        return self.out + '_allframes.tab'

    def translatednaframesfiltering(self, outputranslatednaframes):
        self.df1 = pd.read_csv(outputranslatednaframes, sep="\t", header=None, names=['seqid', 'nseq'])
        self.df2 = self.df1[self.df1['nseq'].str.contains('\*') == False]
        self.df2.to_csv(self.out + '_allframesfiltered.tab', sep="\t",
                        header=None, index=False)
        return self.out + '_allframesfiltered.tab'

    def blastpfilterinf(self, outputnohashblastp):
        self.df1 = pd.read_csv(outputnohashblastp, sep="\t", header=None,
                               names=['chrstartend', 'proteinID', 'percmatch', 'length', 'mismatch', 'gapopen',
                                      'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sseq'])
        # find the best match
        self.df2 = self.df1.groupby(['chrstartend'], as_index=False)[['percmatch', ]].max()
        # primary key for merging
        self.df2['index'] = self.df2['chrstartend'].map(str) + '_' + self.df2['percmatch'].map(str)
        self.df1['index'] = self.df1['chrstartend'].map(str) + '_' + self.df1['percmatch'].map(str)
        #
        self.df3 = pd.merge(self.df1, self.df2, on="index")
        self.df3 = self.df3.loc[self.df3['mismatch'] < 2].reset_index(drop=True)
        self.df3['chr'] = self.df3['chrstartend_x'].str.split(':', 1).str[0]
        self.df3['temp'] = self.df3['chrstartend_x'].str.split(':', 1).str[1]
        self.df3['start'] = self.df3['temp'].str.split('-', 1).str[0]
        self.df3['temp_end'] = self.df3['temp'].str.split('-', 1).str[1]
        self.df3['end'] = self.df3['temp_end'].str.split('.', 1).str[0]
        self.df3['frame'] = self.df3['temp_end'].str.split('.', 1).str[1]
        # self.df3['dupli'] = self.df3.duplicated('index')
        self.df3[['chr', 'start', 'end', 'proteinID', 'frame', 'sstart', 'send', 'sseq']].to_csv(
            self.out + '_blastpnohashfiltered.tab', header=None, sep='\t', index=False)
        return self.out + '_blastpnohashfiltered.tab'

    def outputparsing(self, outputfromaddsequence, outputfromblastpfilterinf):
        self.df1 = pd.read_csv(outputfromaddsequence, sep="\t", header=None,
                               names=['chr', 'clonestart', 'cloneend', 'clonelength', 'start', 'end', 'geneid',
                                      'strand', 'genename', 'description', 'nseq'])
        self.df2 = pd.read_csv(outputfromblastpfilterinf, sep="\t", header=None,
                               names=['chr', 'clonestart', 'cloneend', 'proteinID', 'frame', 'sstart', 'send', 'sseq'])
        self.df3 = pd.merge(self.df1, self.df2, on='clonestart')
        self.df3[['chr_x', 'clonestart', 'cloneend_x', 'clonelength', 'start', 'end', 'geneid', 'strand', 'genename',
                  'description', 'nseq', 'proteinID', 'frame', 'sstart', 'send', 'sseq']].to_csv(
            self.out + '_domaindetection_step1.tab', header=None, sep='\t', index=False)
        return self.out + '_domaindetection_step1_clones_sequence.tab'

    def cleantempfile(self, filedict):
        try:
            for c in filedict.items():
                os.remove(c[1])
        except OSError:
            pass
        else:
            return
