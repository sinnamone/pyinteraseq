import urlparse
from collections import namedtuple
import gzip
import pandas as pd

gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)


def parsegffattributes(attributestring):
    """Parse the GFF3 attribute column and return a dict"""  
    if attributestring == ".":
        return {}
    ret = {}
    for attribute in attributestring.split(";"):
        key, value = attribute.split("=")
        ret[urlparse.unquote(key)] = urlparse.unquote(value)
    return ret


def parsegff3(filename):
    """
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.

    Supports transparent gzip decompression.
    """
    # Parse with transparent decompression
    openfunc = gzip.open if filename.endswith(".gz") else open
    with openfunc(filename) as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            # If this fails, the file format is not standard-compatible
            assert len(parts) == len(gffInfoFields)
            # Normalize data
            normalizedinfo = {
                "seqid": None if parts[0] == "." else urlparse.unquote(parts[0]),
                "source": None if parts[1] == "." else urlparse.unquote(parts[1]),
                "type": None if parts[2] == "." else urlparse.unquote(parts[2]),
                "start": None if parts[3] == "." else int(parts[3]),
                "end": None if parts[4] == "." else int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": None if parts[6] == "." else urlparse.unquote(parts[6]),
                "phase": None if parts[7] == "." else urlparse.unquote(parts[7]),
                "attributes": parsegffattributes(parts[8])
            }
            # Alternatively, you can emit the dictionary here, if you need mutability:
            #    yield normalizedinfo
            yield GFFRecord(**normalizedinfo)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--gff', action="store", dest="gff", default=None,
                        help='Annotation File(.gff) Query at https://www.ncbi.nlm.nih.gov/genome/')
    parser.add_argument("--outputfolder", action="store", dest="outputfolder", default=None, help='Folder output')
    parser.add_argument("--outputid", action="store", dest="outputid", default=None, help='Folder Id')
    args = parser.parse_args()
    if args.outputfolder.endswith('/') is False:
        args.outputfolder = args.outputfolder + '/'
    # Execute the parser
    templist = []
    # iterate in file gff converted in dictionary
    for record in parsegff3(args.gff):
        templist.append(record)
    # copy dictionary in dataframe
    df = pd.DataFrame(templist)
    # split dataframe in two dataframe, df1a filter gene df1b filter CDS
    df1a = df.loc[df['type'] == 'gene'].reset_index(drop=True)
    # split columns attributes in 7 columns
    df2a = df1a['attributes'].apply(pd.Series)
    # index of new db became start
    df2a['id'] = df1a['start']
    # subset dataframe to get gene id that is labeled as locus_tag
    df3a = df2a[['id', 'locus_tag']]
    # merge
    df4a = pd.merge(df1a, df3a, left_on='start', right_on='id')
    # split dataframe in two dataframe, df1a filter gene df1b filter CDS
    df1b = df.loc[df['type'] == 'CDS'].reset_index(drop=True)
    # split columns attributes in 7 columns
    df2b = df1b['attributes'].apply(pd.Series)
    # index of new db became start
    df2b['id'] = df1b['start']
    # subset dataframe to get description id that is labeled as product
    df3b = df2b[['id', 'product']]
    # merge
    df4b = pd.merge(df1b, df3b, left_on='start', right_on='id')
    # merge all information
    df5 = pd.merge(df4b, df4a, on='start')
    df5['score'] = "."
    # save in tsv
    df5[['seqid_x', 'start', 'end_x', 'locus_tag', 'score', 'strand_x', 'product']].to_csv(
        args.outputfolder + args.outputid + '.bed', sep="\t", header=None, index=False)
