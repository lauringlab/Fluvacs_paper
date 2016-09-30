import argparse
import glob
import pandas as pd
import os

parser = argparse.ArgumentParser(description='This script is designed to parse the concatenated consensus sequence from the deepSNV output into 1) a fastq file that has a sequence for each segment.')


parser.add_argument('fasta', metavar='dir', nargs='+',
                    help='The fasta file that is to be parsed.')

parser.add_argument('key', metavar='key', nargs='+',
                    help='This is the csv that contains the concatenated positons of the segments (base 1)')

parser.add_argument('out_fa', metavar='out_dir', nargs='+',
                    help='The output file that is to be made')


args = parser.parse_args()

fa=args.fasta[0]
key=args.key[0]
out_file=args.out_fa[0]

def get_seg(string,segment,seg_df):
    "This splits a sting by position into the 8 genomic segments"
    seg=seg_df['chr'][seg_df['chr']==segment].index[0] # get the index of the segment
    start=seg_df['start'][seg]-1 # the positions are 1 based but we need them 0 based
    stop=seg_df['stop'][seg]-1
    seq=string[start:stop]
    return(seq)
def print_seg(string,segment,seg_df,outfile):
    "This writes the segment to a file"
    seq=get_seg(string,segment,seg_df)
    with open(outfile,'a') as out: # a so that we append each time
        out.write(">"+segment+"\n"+seq+"\n")


seg_df=pd.read_csv(key, sep=',',comment='#')

segs=seg_df['chr']


print 'working with ' + fa
with open(fa,'r') as f:
    line=f.readlines()[0] # There is only 1 line
for chr in segs:
    #print chr
    print_seg(line,chr,seg_df,out_file)
