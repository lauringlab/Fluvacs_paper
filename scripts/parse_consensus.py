#!/Users/jt/.virtualenvs/sci-py3/bin/python
import argparse
import glob
import pandas as pd
import os

parser = argparse.ArgumentParser(description='This script is designed to parse the concatenated consensus sequence from the deepSNV output into 1) a fastq file that has a sequence for each segment or 2) one file that has the sequence of a given segment from each sample.\n By default the script takes in a directory and works with all fasta files.')


parser.add_argument('dir', metavar='dir', nargs='+',
                    help='The directory that contains the fasta files')

parser.add_argument('key', metavar='key', nargs='+',
                    help='This is the csv that contains the concatenated positons of the segments (base 1)')

parser.add_argument('out_dir', metavar='out_dir', nargs='+',
                    help='The output directory. If it doesn\'t exist it will be made')

parser.add_argument('--seg', dest='segment', action='store',
                    default= 'none', help='The segment to be taken from each file' )

args = parser.parse_args()

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
def print_seg_file(string,segment,seg_df,outfile,name):
    "This writes the segment and file name to a file "
    seq=get_seg(string,segment,seg_df)
    with open(outfile,'a') as out: # a so that we append each time
        out.write(">"+name+"\n"+seq+"\n")

seg_df=pd.read_csv(args.key[0], sep=',',comment='#')

## If the output dir does not exist make it
if not os.path.exists(args.out_dir[0]):
    os.makedirs(args.out_dir[0])


segs=seg_df['chr']
if args.segment=='none':
    for fa in glob.glob(args.dir[0]+"/*.fa"):
        print('working with ' + fa)
        with open(fa,'r') as f:
            line=f.readlines()[0] # There is only 1 line
        for chr in segs:
            #print(chr)
            print_seg(line,chr,seg_df,args.out_dir[0]+"/"+os.path.basename(fa))
if args.segment!='none':
    for fa in glob.glob(args.dir[0]+'/*.fa'):
        print('working with ' + fa)
        with open(fa,'r') as f:
            line=f.readlines()[0] # There is only 1 line
        print_seg_file(line,args.segment,seg_df,args.out_dir[0]+'/'+args.segment+'.fa',os.path.basename(fa))
