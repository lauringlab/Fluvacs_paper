#import sys
import os
import gzip
import argparse
import shutil
import glob
parser = argparse.ArgumentParser(description='This program takes hiseq fastq files and renames them as sample.read_direction.#.fastq and keeps a log of the change. It renames the files according the Library name in the runifo csv. This relies heavily on the change_name scripts in Variant_pipeline. -JT McCrone 2016')
parser.add_argument('-s',action='store',dest='s',help='The sorce directory containing the original fastq files')
parser.add_argument('-f',action='store',dest='f',help='The final directory that will hold the renamed fastq files')
parser.add_argument('-k',action='store',dest='key',help='The sra runinfo csv that serves as the key for renaming')
parser.add_argument('-run',action='store_true',dest='test',default=False,help='Boolean switch to run program, without this the program runs in test mode: the log is made but no files are renamed')

args=parser.parse_args()
s=args.s
f=args.f
key=args.key
test=args.test

#make the output dir if it doesn't exist
if not os.path.exists(f):
    os.makedirs(f)

# add slash to final location if it was forgotten
if f[-1] != '/':
    f=f+'/'
#print f

if s[-1] != '/':
    s=s+'/'
#print s


# input argument is the sample sheet
junk_names = [] # The names given by the sequencing core sampleid_index
new_names = [] # the names we want, which are in the same row as the data taht is used to make the bad names. 

# so junk_names[i] should be the junk name for new_names[i] as they both come from the ith column.
# add the bad names from the sample sheet to a list. This uses columns of the sheet to construct the names according to the capricous nature of the sequencing core.Although they may have finally settled down a bit.
#print(key)
names =  open(key,"r")
next(names) # skip the header of the csv
for line in names:
    line = line.strip()
    line = line.split(',')
    junk_names.append(line[0])    
    new=line[1]   
    new_names.append(new)
   # print(line)
names.close()

#print(new_names)
#print(junk_names)
if test==False:
    print "running in test mode add option -run to run"
    
outfile = open(f+'renaming_log.txt','w')


# Now to search through the fastq files or fastq.gz and copy them to the new makes
for filename in glob.glob(s + "*.fastq"):
    #print(filename)
    path=os.path.abspath(filename) # getteh full name and path
    name=filename.split("/")[-1]
    bad_name = name.split("_")[0] # the bad name is the first bit
    #print(bad_name)
    lane_junk = name.split("_")[1]
    read_number=lane_junk.split(".fastq") # info about the read 
    fastq_number=str(1) #always 1 from the sra
    read_number=str(read_number[0]) # read direction
    if bad_name in junk_names: # the new name and the bad name have the same index in their respective lists
        name_index = junk_names.index(bad_name)
        better_name= new_names[name_index]
        perfect_name= better_name+"."+read_number+"."+fastq_number+".fastq"
        # Write file to new name
        print("Moving "+ path + " to "+f+perfect_name)
        outfile.write(path + "\t Moving to \t" + f+perfect_name + "\n") # log the move in a renaming_log.txt
        if test==True:
            shutil.move(path,f+perfect_name)
            
#print(junk_names)
#print(new_names)


outfile.close()

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
# Make sure the files were moved and delete the log if they weren't

if file_len(f+'renaming_log.txt')<1:
    os.remove(outfile)
    print("Renaming failed - less than 1 line in the log")
