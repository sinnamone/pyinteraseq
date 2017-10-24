import subprocess
from output_message import *
from Bio import SeqIO
import os
import pandas as pd
import sys
from pyinteraseq_inputcheck import InputCheck



class BlastNlucleotide(InputCheck):

    def __init__(self, optparseinstance):
        InputCheck.__init__(self, optparseinstance)
        self.out_lines = []
        self.temp_line = ''
        self.df1 = None
        self.seqix = 1
        self.id = '1'
        self.id = [(int(x) - 1) for x in self.id]
        self.header = None
        self.dbid = os.path.basename(self.fastasequence.split('/')[-1])
        self.blastnformat7 = '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq'
        self.blastnformat6 = '6 sseqid sstart send qseqid score sstrand'

    def fastq2fasta(self, fastq, nameid):
        if nameid == "forward":
            SeqIO.convert(fastq, 'fastq', self.out + '1.fasta', 'fasta')
            return self.out + '1.fasta'
        elif nameid == "reverse":
            SeqIO.convert(fastq, 'fastq', self.out + '1.fasta', 'fasta')
            return self.out + '1.fasta'

    def fasta2tab(self, fasta, nameid):
        if nameid == "forward":
            with open(fasta, 'r') as fp:
                for line in fp:
                    if line.startswith('>'):
                        self.out_lines.append(self.temp_line)
                        self.temp_line = line.strip() + '\t'
                    else:
                        self.temp_line += line.strip()
            with open(self.out + '1.tab', 'w') as fp_out:
                fp_out.write('\n'.join(self.out_lines))
            with open(self.out + '1.tab', 'r+') as f:  # open in read / write mode
                f.readline()  # read the first line and throw it out
                data = f.read()  # read the rest
                f.seek(0)  # set the cursor to the top of the file
                f.write(data)  # write the data back
                f.truncate()  # set the file size to the current size
            fp_out.close()
            f.close()
            return self.out + '1.tab'
        elif nameid == "reverse":
            with open(fasta, 'r') as fp:
                for line in fp:
                    if line.startswith('>'):
                        self.out_lines.append(self.temp_line)
                        self.temp_line = line.strip() + '\t'
                    else:
                        self.temp_line += line.strip()
            with open(self.out + '2.tab', 'w') as fp_out:
                fp_out.write('\n'.join(self.out_lines))
            with open(self.out + '2.tab', 'r+') as f:  # open in read / write mode
                f.readline()  # read the first line and throw it out
                data = f.read()  # read the rest
                f.seek(0)  # set the cursor to the top of the file
                f.write(data)  # write the data back
                f.truncate()  # set the file size to the current size
            fp_out.close()
            f.close()
            return self.out + '2.tab'
        elif nameid == "domains":
            with open(fasta, 'r') as fp:
                for line in fp:
                    if line.startswith('>'):
                        self.out_lines.append(self.temp_line)
                        self.temp_line = line.strip() + '\t'
                    else:
                        self.temp_line += line.strip()
            with open(self.out + '_domain.tab', 'w') as fp_out:
                fp_out.write('\n'.join(self.out_lines))
            with open(self.out + '_domain.tab', 'r+') as f:  # open in read / write mode
                f.readline()  # read the first line and throw it out
                data = f.read()  # read the rest
                f.seek(0)  # set the cursor to the top of the file
                f.write(data)  # write the data back
                f.truncate()  # set the file size to the current size
            fp_out.close()
            f.close()
            return self.out + '_domains.tab'

    def seqrename(self, tabular, readirection):
        print readirection
        if readirection == "forward":
            self.df1 = pd.read_csv(tabular, header=None, sep='\t')
            self.df1['seq_id'] = self.df1.apply(lambda x: "seq1:1:" + str(x.name), axis=1)
            self.df1[['seq_id', 1]].to_csv(self.out + '_1_newid.tab', header=None, sep='\t', index=False)
            return self.out + '_1_newid.tab'
        else:
            self.df1 = pd.read_csv(tabular, header=None, sep='\t')
            self.df1['seq_id'] = self.df1.apply(lambda x: "seq2:2:" + str(x.name), axis=1)
            self.df1[['seq_id', 1]].to_csv(self.out + '_2_newid.tab', header=None, sep='\t', index=False)
            return self.out + '_2_newid.tab'

    def tab2fasta(self, tabular, readirection):
        print readirection
        with open(tabular, 'r') as f:
            if readirection == "forward":
                with open(self.out + '_1_newid.fasta', 'w') as f_out:
                    for line in f:
                        line = line.strip().split('\t')
                        self.header = '>' + '_'.join([line[i] for i in self.id])
                        f_out.write(self.header + '\n')
                        f_out.write(line[self.seqix] + '\n')
                f_out.close()
                return self.out + '_1_newid.fasta'
            elif readirection == "reverse":
                with open(self.out + '_2_newid.fasta', 'w') as f_out:
                    for line in f:
                        line = line.strip().split('\t')
                        self.header = '>' + '_'.join([line[i] for i in self.id])
                        f_out.write(self.header + '\n')
                        f_out.write(line[self.seqix] + '\n')
                f_out.close()
                return self.out + '_2_newid.fasta'
            elif readirection == "single":
                with open(self.out + '_newid.fasta', 'w') as f_out:
                    for line in f:
                        line = line.strip().split('\t')
                        self.header = '>' + '_'.join([line[i] for i in self.id])
                        f_out.write(self.header + '\n')
                        f_out.write(line[self.seqix] + '\n')
                f_out.close()
                return self.out + '_newid.fasta'

    def concatenateforrev(self, readlist):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        with open(self.out + '_con.fasta', 'w') as outfile:
            for fname in readlist:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
        outfile.close()
        # self.filelog.write(msg56 + self.fastqcount(self.out + '_con.fasta',
        #                                            self.readforwardtype))
        return self.out + '_con.fasta'

    def makeblastnucl(self, fasta):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            if os.path.isfile(self.outputfolder + self.dbid + '.db') is False:
                subprocess.check_call(
                    ['makeblastdb',
                     '-in', fasta,
                     '-dbtype',
                     'nucl',
                     '-out', self.outputfolder + self.dbid + 'blastn.db'],
                    stderr=self.filelog, stdout=self.filelog)
            else:
                self.filelog.write(msg57)
        except subprocess.CalledProcessError:
            self.filelog.write(msg58)
            sys.exit(0)
        else:
            self.filelog.write(msg59)
            return self.outputfolder + self.dbid + 'blastn.db'

    def blastn(self, fastainpu, fasta):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            subprocess.check_call(['blastn', '-out', self.out + '_blastn.tab',
                                   '-outfmt',
                                   self.blastnformat7,
                                   '-query', fastainpu,
                                   '-db', self.makeblastnucl(fasta),
                                   '-evalue', '0.001',
                                   '-num_threads', self.thread],
                                  stderr=self.filelog, stdout=self.filelog)
        except subprocess.CalledProcessError:
            self.filelog.write(msg60)
            sys.exit(0)
        else:
            self.filelog.write(msg61)
            return self.out + '_blastn.tab'

    def blastnclones(self, fastainpu, fasta):
        self.filelog = open(self.outputfolder + self.outputid + ".log", "a")
        try:
            subprocess.check_call(['blastn', '-out', self.out + '_blastnclones.tab',
                                   '-outfmt',
                                   self.blastnformat6,
                                   '-query', fastainpu,
                                   '-db', self.makeblastnucl(fasta),
                                   '-evalue', '0.001',
                                   '-num_threads', self.thread],
                                  stderr=self.filelog, stdout=self.filelog)
        except subprocess.CalledProcessError:
            self.filelog.write(msg60)
            sys.exit(0)
        else:
            self.filelog.write(msg61)
            return self.out + '_blastnclones.tab'

    def fasta2tabular(self, imp, prefix):
        SeqIO.convert(imp, 'fasta', self.out + prefix + '.tab', 'tab')
        return self.out + prefix + '.tab'





