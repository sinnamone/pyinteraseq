import urlparse
from collections import namedtuple
import gzip
import pandas as pd
import os
import traceback
import optparse
from Bio import SeqIO



class InputParsing(object):

    def __init__(self, optparseinstance):
        self.inputistance = optparseinstance
        # read input flag
        #self.fastasequence = self.inputistance.fastasequence
        self.gff = self.inputistance.gff
        self.outputfolder = self.inputistance.outputfolder
        self.outputid = self.inputistance.outputid
        #
        self.gffInfoFields = ['Location', 'Strand', 'Length', 'PID', 'Gene', 'Synonym', 'Code', 'COG', 'Product']
        self.GFFRecord = namedtuple("GFFRecord", self.gffInfoFields)
        # name of fasta file without path
        #self.namefilefasta = os.path.basename(self.fastasequence.split('/')[-1])
        # path were fasta file is located
        #self.PathFasta = os.path.dirname(self.fastasequence)
        #self.chromosomename = self.namefilefasta.split('.')[0]
        self.genome = None
        self.ref = None
        self.index = None
        self.numchr = None
        self.df1 = None
        self.dir = os.listdir(self.outputfolder)
        self.id = '1'
        self.id = [(int(x) - 1) for x in self.id]
        self.header = None
        self.seqix = 1
        # for self.seq_record in SeqIO.parse(self.fastasequence, "fasta"):
        #     # save the name of chromosome. Format string is composed by >gi|XXXXXX|ref|NC_XXXXX.1|description
        #     # self.chromosomename = self.seq_record.id.split("|")[-2]
        #     self.chromosomelength = len(self.seq_record)
        # print self.chromosomename
        #
        if self.outputfolder is not None:
            if self.outputfolder.endswith('/') is True:
                self.outputfolder = self.inputistance.outputfolder
                self.out = self.outputfolder + self.outputid
            else:
                self.outputfolder = self.inputistance.outputfolder + '/'
                self.out = self.outputfolder + self.outputid

    def openlog(self):
        return open(self.outputfolder + self.outputid + "_inputparsing.log", "a")

    def genomefilewrite(self):
        """
        Create file .genome
        :return:
        """
        self.filelog = self.openlog()
        try:
            with open(self.outputfolder + self.chromosomename + '.genome', "wb") as self.genome:
                print self.genome
                self.genome.write(self.chromosomename + '\t' + str(self.chromosomelength))
        except OSError:
            self.filelog.write('msg45')
            sys.exit(0)
        else:
            self.filelog.write('msg46')
            return self.outputfolder + self.chromosomename + '.genome'

    def ref2tabular(self):
        """
        Function for conversion Reference Fasta file in tabular format
        :return:
        """
        self.filelog = self.openlog()
        fileout = self.outputfolder + self.chromosomename + '.tab'
        try:
            SeqIO.convert(self.fastasequence, 'fasta', fileout, 'tab')
        except traceback:
            self.filelog.write(traceback.format_exc())
            self.filelog.write('msg8')
            sys.exit(1)
        else:
            self.filelog.write('msg7')
            return fileout

    # def chrcount(self):
    #     """
    #     Rename chr id of reference fasta
    #     :return:
    #     """
    #     self.filelog = self.openlog()
    #     try:
    #         # output variable
    #         outid = self.outputfolder + self.chromosomename + '.tab'
    #         # count number of chr of fasta file
    #         self.numchr = sum(1 for line in open(self.ref2tabular()))
    #         # create dataframe converting fasta to tabular
    #         self.df1 = pd.read_csv(self.ref2tabular(), header=None, sep='\t')
    #         pippo self.df1[0].to_dict(orient="index")
    #         # change id chr
    #         self.df1[0] = self.df1.apply(lambda x: self.chromosomename + '_' + str(x.name + 1), axis=1)
    #         # output write
    #         self.df1.to_csv(outid, header=None, sep='\t', index=False)
    #     except traceback:
    #         self.filelog.write(traceback.format_exc())
    #         self.filelog.write('msg114')
    #         sys.exit(1)
    #     else:
    #         self.filelog.write('msg115')
    #         return outid

    def fastareference(self):
        self.filelog = self.openlog()
        try:
            with open(self.outputfolder + self.chromosomename + '.fasta', 'w') as self.ref:
                print self.ref
                for self.seq_record in SeqIO.parse(self.fastasequence, "fasta"):
                    self.ref.write('>' + str(self.chromosomename) + '\n' + str(self.seq_record.seq))
        except OSError:
            self.filelog.write(msg43)
            sys.exit(1)
        else:
            for item in self.dir:
                if item.endswith(".fai"):
                    os.remove(self.outputfolder + item)
            self.filelog.write(msg44)
            self.index = pysam.Fastafile(self.outputfolder + self.chromosomename + '.fasta')
            print self.index
            return self.outputfolder + self.chromosomename + '.fasta'

    # def tab2fasta(self):
    #     """
    #     Convert tabular fasta file (2 columns) first ID second nucleotide sequence in fasta file format
    #     :return:
    #     """
    #     self.filelog = self.openlog()
    #     try:
    #         with open(self.out + self.chromosomename + '.tab', 'r') as f:
    #             with open(self.out + 'reformat.fasta', 'w') as f_out:
    #                 for line in f:
    #                     line = line.strip().split('\t')
    #                     self.header = '>' + '_'.join([line[i] for i in self.id])
    #                     f_out.write(self.header + '\n')
    #                     f_out.write(line[self.seqix] + '\n')
    #             f_out.close()
    #     except traceback:
    #         self.filelog.write(traceback.format_exc())
    #         self.filelog.write('msg100')
    #         sys.exit(1)
    #     else:
    #         self.filelog.write('msg101')
    #         return self.out + 'reformat.fasta'

    def parsegffattributes(self, attributestring):
        """Parse the GFF3 attribute column and return a dict"""
        if attributestring == ".":
            return {}
        ret = {}
        for attribute in attributestring.split(";"):
            key, value = attribute.split("=")
            ret[urlparse.unquote(key)] = urlparse.unquote(value)
        return ret

    def parsegff3(self, filename):
        """
        A minimalistic GFF3 format parser.
        Yields objects that contain info about a single GFF3 feature.

        Supports transparent gzip decompression.
        """
        # Parse with transparent decompression
        openfunc = gzip.open if filename.endswith(".gz") else open
        with openfunc(filename) as infile:
            for _ in xrange(3):
                next(infile)
            for line in infile:
                parts = line.strip().split("\t")
                # If this fails, the file format is not standard-compatible
                assert len(parts) == len(self.gffInfoFields)
                # Normalize data
                normalizedinfo = {
                    "Location": None if parts[0] == "." else urlparse.unquote(parts[0]),
                    "Strand": None if parts[1] == "." else urlparse.unquote(parts[1]),
                    "Length": None if parts[2] == "." else int(parts[2]),
                    "PID": None if parts[3] == "." else int(parts[3]),
                    "Gene": None if parts[4] == "." else urlparse.unquote(parts[4]),
                    "Synonym": None if parts[5] == "." else urlparse.unquote(parts[5]),
                    "Code": None if parts[6] == "." else urlparse.unquote(parts[6]),
                    "COG": None if parts[7] == "." else urlparse.unquote(parts[7]),
                    "Product": None if parts[8] == "." else urlparse.unquote(parts[8])
                }
                # Alternatively, you can emit the dictionary here, if you need mutability:
                #    yield normalizedinfo
                yield self.GFFRecord(**normalizedinfo)







if __name__ == "__main__":
    parser = optparse.OptionParser(usage='python %prog', version='1.0', )
    # parser.add_option('--fastasequence', action="store", dest="fastasequence", default=None,
    #                   help='Read dataset input forward')
    parser.add_option('--gff', action="store", dest="gff", default=None,
                      help='Read dataset input reverse')
    parser.add_option("--outputfolder", action="store", dest="outputfolder", default=None, help='Folder output')
    parser.add_option("--outputid", action="store", dest="outputid", default=None, help='Folder Id')
    options, args = parser.parse_args()
    DictInfo = dict()
    ParserClass = InputParsing(optparseinstance=options)

    # Execute the parser
    templist = []  # iterate in file gff converted in dictionary
    for record in ParserClass.parsegff3(options.gff):
        templist.append(record)
    # copy dictionary in dataframe
    dfAnno = pd.DataFrame(templist)
    aux = dfAnno['Location'].apply(lambda x: x.split('..'))
    dfAnno['Start'] = aux.apply(lambda x: x[0])
    dfAnno['End'] = aux.apply(lambda x: x[1])
    dfAnno['Chr'] = os.path.basename(options.gff.split('/')[-1]).split('.')[0]
    dfAnno.loc[dfAnno['Strand'] == '+', ['Start']] = dfAnno['Start'].astype(int).sub(1)
    dfAnno[['Chr', 'Start', 'End', 'Synonym', 'Code', 'Strand', 'Product']].to_csv(
        options.outputfolder + options.outputid + '.bed', header=None, sep='\t', index=False)
