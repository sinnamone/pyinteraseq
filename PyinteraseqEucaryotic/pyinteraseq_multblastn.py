from itertools import izip_longest
import subprocess
from multiprocessing import Pool
import os
import optparse
import sys
from output_message import *

parser = optparse.OptionParser(usage='python %prog pyinteraseq_multiblastn.py', version='1.0',)
parser.add_option('--referencefasta', action="store", dest="referencefasta", default=None,
                  help='Genome sequence fasta')
parser.add_option('--multifastasequence', action="store", dest="multifastasequence", default=None,
                  help='Trimmed dataset in multifasta format.')
parser.add_option('--dbname', action="store", dest="dbname", default=None,
                  help='Name given to makeblastdb database.')
parser.add_option('--chunks', action="store", dest="chunks", default=50000,
                  help='Output target derived from pyinteraseq_mapping.py',)

query_opts = optparse.OptionGroup(
    parser, 'Output Options',
    'Options for the output destionation and name.',
    )
query_opts.add_option('--outputfolder', action="store", dest="outputfolder", default=None,
                      help='Output folder.')
query_opts.add_option('--outputid', action="store", dest="outputid", default=None,
                      help='Output ID.')
query_opts.add_option('--outformat', action="store", dest="outformat", default=None,
                      help='Output format for blastn.')
query_opts.add_option('--suffix', action="store", dest="suffix", default='_blastn.txt',
                      help='Suffix to add after merging files.')
parser.add_option_group(query_opts)

reference_opts = optparse.OptionGroup(
    parser, 'Advanced Options',
    'Options for advanced analysis.',
    )
reference_opts.add_option('--thread', action="store", dest="thread", default='2',
                          help='Number of thread.')
parser.add_option_group(reference_opts)

options, args = parser.parse_args()


def grouper(ni, iterable, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks
    :param ni:
    :param iterable:
    :param fillvalue:
    :return:
    """
    args = [iter(iterable)] * ni
    return izip_longest(fillvalue=fillvalue, *args)


def splittingfiles(wholefasta, chunks, outputpath, idtemp):
    """
    Split multifasta sequence files
    :param wholefasta: Input Fasta file
    :param chunks: line number of temp fasta file
    :param outputpath: Output path
    :param idtemp: Prefix to give at temp file
    :return: number of file generated
    """
    with open(wholefasta) as f:
        count = 0
        for i, g in enumerate(grouper(chunks, f, fillvalue=''), 1):
            count += 1
            with open(outputpath + '{0}_temp_{1}.fasta'.format(idtemp, i * chunks), 'w') as fout:
                fout.writelines(g)
    return count


def makeblastdb(outputfolder, dbname, fasta):
    """
    Run makeblastdb
    :param outputfolder: Output folder for new db
    :param dbname: Name of Database
    :param fasta: Reference fasta file
    :return: Database name
    """
    if os.path.isfile(outputfolder + dbname) is False:
        fnull = open(os.devnull, 'w')
        try:
            subprocess.check_call(
                ['/opt/ncbi-blast-2.7.1+/bin/makeblastdb',
                 '-in', fasta,
                 '-dbtype',
                 'nucl',
                 '-out', options.outputfolder + options.dbname],
                stdout=fnull, stderr=fnull)
        except subprocess.CalledProcessError:
            sys.exit(1)
        else:
            return outputfolder + dbname


def blastn(outputname, fastainpu, dbname, outputformat):
    """
    Run blastn
    :param outputname: Prefix given to temp blastn file
    :param fastainpu: Temp fasta input
    :param dbname: Database name
    :param outputformat: outpformat blast
    :return:
    """
    fnull = open(os.devnull, 'w')
    return (subprocess.check_call(['/opt/ncbi-blast-2.7.1+/bin/blastn', '-out', outputname,
                                   '-outfmt',
                                   outputformat,
                                   '-query', fastainpu,
                                   '-db', dbname,
                                   '-evalue', '0.001'],
                                  stdout=fnull, stderr=fnull))


def blastmultiprocess(maxchunks, outmerge, chunks, thread, databasename, outputformat):
    """
    Controller to run multithread blastn
    :param maxchunks: Max chunk
    :param outmerge: Output Path + Output Prefix
    :param chunks: Number of line for temp file
    :param thread: Number of processor
    :param databasename: Database name
    :param outputformat: blstn output format (6=tab,7= hash)
    :return: List with two list, first temp blastn file, second tempo fasta file
    """
    listempfiles = []
    listfastatemp = []
    pool = Pool(processes=int(thread))
    for i in range(chunks, maxchunks, chunks):
        pool.apply_async(blastn, args=(outmerge + '_{0}.txt'.format(i),
                                       outmerge + '_temp_{0}.fasta'.format(i),
                                       databasename, outputformat),)
        listempfiles.append(outmerge + '_{0}.txt'.format(i))
        listfastatemp.append(outmerge + '_temp_{0}.fasta'.format(i))
    pool.apply_async(blastn, args=(outmerge + '_{0}.txt'.format(maxchunks),
                                   outmerge + '_temp_{0}.fasta'.format(maxchunks),
                                   databasename, outputformat),)
    listempfiles.append(outmerge + '_{0}.txt'.format(maxchunks))
    listfastatemp.append(outmerge + '_temp_{0}.fasta'.format(maxchunks))
    pool.close()
    pool.join()
    return listempfiles, listfastatemp


def mergtempfiles(listtemporanyfiles):
    """
    Merge temporany files
    :param listtemporanyfiles: list created by blastmultiprocess
    :return: absolute path + merged file
    """
    with open(options.outputfolder + options.outputid, 'w') as outfile:
        for fname in listtemporanyfiles:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
    return listtemporanyfiles


def cleantempfiles(listtemporanyfiles):
    """
    Remove temporary file
    :param listtemporanyfiles:
    :return:
    """
    for fname in listtemporanyfiles:
        os.remove(fname)
    return listtemporanyfiles


if __name__ == '__main__':
    # Check input path
    if options.outputfolder is not None:
        if options.outputfolder.endswith('/') is True:
            outp = options.outputfolder + options.outputid
        else:
            outp = options.outputfolder + '/' + options.outputid
    # Makeblastdb
    databname = makeblastdb(outputfolder=options.outputfolder, dbname=options.dbname, fasta=options.referencefasta)
    # Split input multifastafile
    numfilesgenerated = (splittingfiles(wholefasta=options.multifastasequence, chunks=options.chunks,
                                        outputpath=options.outputfolder,
                                        idtemp=options.outputid))*options.chunks
    # blastn on multi thread
    blastntemp = blastmultiprocess(maxchunks=numfilesgenerated, outmerge=outp, chunks=options.chunks,
                                   thread=options.thread, databasename=databname,
                                   outputformat=options.outformat)
    # Merge temp file
    mergtempfiles(listtemporanyfiles=blastntemp[0])
    # Remove temp file
    for i in range(2):
        cleantempfiles(listtemporanyfiles=blastntemp[i])
