import optparse
from Bio import SeqIO
import os
import pandas as pd
import sys

parser = optparse.OptionParser(usage='python %prog -i [ input file with path ] -o [ output ID ] -p [ output path ] -t [ selection between \'fastq2fasta\',\'fasta2tab\',\'seqrename\',\'tab2fasta\' ] -r [ 1 or 2 ]',version='1.0',)
parser.add_option('-i', action="store", dest="input",help='input fasta,fastq,tabular file included path')
parser.add_option('-o', action="store", dest="output",help='Output Id')
parser.add_option('-p', action="store", dest="pathout",help='Output Pathway')
parser.add_option('-t', type='choice', choices=['fastq2fasta','fasta2tab','seqrename','tab2fasta'],help='Insert fastq2fasta for conversion fastq to fasta; Insert fasta2tab per conversion to Fasta to Tabular; Insert seqrename for change read name with generic seq+line num, Insert tab2fasta to convert tabular file (2 columns, firt with reads name without > and second with sequence.)')
parser.add_option('-r', type='choice', choices=['1','2'], default= '1',help='1=forward;2=reverse')
options, args = parser.parse_args()

class confastq2fasta(object):
	def __init__(self,input,pathout,output,read):
		self.input = input
		self.output = output
		self.read = read 
		if pathout.endswith('/') == True:
			self.path_out = pathout
		else:
			self.path_out = pathout+'/'

            
	def fastq2fasta(self):
		SeqIO.convert(self.input,'fastq',self.path_out+self.output+'_'+self.read+'.fasta','fasta')
		
class confasta2tab(confastq2fasta):
	def __init__(self,input,pathout,output,read):
		confastq2fasta.__init__(self, input,pathout,output,read)
		self.out_lines = []
		self.temp_line = ''

		
	def fasta2tab(self):
		with open(self.input,'r') as fp:
			for line in fp:
				if line.startswith('>'):
					self.out_lines.append(self.temp_line)
					self.temp_line = line.strip() + '\t'
				else:
					self.temp_line += line.strip()
 		with open(self.path_out+self.output+'_'+self.read+'.tab', 'w') as fp_out:
			fp_out.write('\n'.join(self.out_lines))
		with open(self.path_out+self.output+'_'+self.read+'.tab', 'r+') as f: #open in read / write mode
			f.readline() #read the first line and throw it out
			data = f.read() #read the rest
			f.seek(0) #set the cursor to the top of the file
			f.write(data) #write the data back
			f.truncate() #set the file size to the current size

class countseq(confasta2tab):
	def __init__(self,input,pathout,output,read):
		confasta2tab.__init__(self, input,pathout,output,read)
		
		
		
	def seqrename(self):
		self.df1=pd.read_csv(self.input,header=None,sep='\t')
		self.df1['seq_id'] =  self.df1.apply(lambda x: "seq"+self.read+"_" + str(x.name), axis=1)
		self.df1[['seq_id',1]].to_csv(self.path_out+self.output+'_'+self.read+'_newid.tab',header=None, sep='\t', index=False)


class contabular2fasta(countseq):
	def __init__(self,input,pathout,output,read):
		countseq.__init__(self,input,pathout,output,read)
		self.seqix = 1
		self.id= '1'
		self.id= [(int(x) - 1) for x in self.id]

	
	def tab2fasta(self):
		with open(self.input,'r') as f:
			with open(self.path_out+self.output+'_'+self.read+'_newid.fasta', 'w') as f_out:
				for line in f:
					line= line.strip().split('\t')
					self.header= '>' + '_'.join([line[i] for i in self.id ])
					f_out.write(self.header+'\n')
					f_out.write(line[self.seqix]+'\n')
			f_out.close()
		f.close()

if __name__ == '__main__':
	if options.t == "fastq2fasta":
		print('Conversion : fastq2fasta')
		a  = confastq2fasta(options.input,options.pathout,options.output,options.r).fastq2fasta()
	elif options.t == "fasta2tab":
		print('Conversion : fasta2tab')
		b  = confasta2tab(options.input,options.pathout,options.output,options.r).fasta2tab()
	elif options.t == "seqrename":
		print('Conversion : seqrename')
		if options.r == "1":
			c  = countseq(options.input,options.pathout,options.output,options.r).seqrename()
		elif options.r == "2":
			c  = countseq(options.input,options.pathout,options.output,options.r).seqrename()
	elif options.t == "tab2fasta":
		print('Conversion : tab2fasta')
		if options.r == "1":
			d  = contabular2fasta(options.input,options.pathout,options.output,options.r).tab2fasta()
		elif options.r == "2":
			d  = contabular2fasta(options.input,options.pathout,options.output,options.r).tab2fasta()
