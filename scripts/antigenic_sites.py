#!/Users/jt/.virtualenvs/sci-py2.7/bin/python

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
import pandas as pd
import os
import tempfile
import sys
import subprocess
import argparse



parser=argparse.ArgumentParser(description="The goal of this work is to determine which HA variants are nonsymonous and whether or any of the nonsymonous variants are in putative antigenic regions. I am only looking at minor variants in this analyis and so nonsymonous variants will be identified as variants that change the amino acid a given position compared to the sample's consensus seqeunce (not the plasmid control).")

parser.add_argument('muscle', metavar='muscle', nargs='+',
                    help='The directory that contains the muscle exicutable')

parser.add_argument('fasta',metavar='fasta',nargs='+',
			help = "The file containing the fasta files of all the samples to be analyzed.")

parser.add_argument('snv',metavar='snv',nargs='+',help = 'The csv containing the snv calls to analyze.')

parser.add_argument('prefix',metavar='prefix',nargs='+',
			help = "The prefix used to save the putative.csv and aa.csv files")

args=parser.parse_args()

muscle_dir=os.path.abspath(args.muscle[0])
fasta=os.path.abspath(args.fasta[0])
snv=os.path.abspath(args.snv[0])
Ids=[]
nucleotide=[]
protein=[]

"""
The first step will be to read in the parsed HA consensus sequences. Here I am reading each sequence translating, and then saving the
sequence names, nucleotide seq object, and protein seq object in separate lists. I'm not inlcuding a stop codon
character in the translation.
"""


for seq_record in SeqIO.parse(fasta, "fasta"):
    Ids.append(seq_record.id.split("_")[0])
    
    seq_record.seq.alphabet=generic_dna
    
    # in sample 427 there is a gap in the consensus near the end of the segment. This sample has low coverage. The consensus reads CA-TTA In every other sample this position is a T. I am making it a T here for translation purposes. It is the only gap in all the HA consesus sequences.
    if seq_record.id.split("_")[0]=='427':
        str_seq=str(seq_record.seq).replace("-","N")
        seq=Seq(str_seq,generic_dna)
        nucleotide.append(seq)
        prot=seq.translate(to_stop=True)
        protein.append(prot)
    else:
        nucleotide.append(seq_record.seq)
        prot=seq_record.seq.translate(to_stop=True)
        protein.append(prot)
"""
Now I will read in all polymorphisms on HA from the 2007-2008 season.
Now to find the codon of each of these positions we will just need to use a little modular algebra.
bp 0-2 are in represent the o aa 3-5 the 1 and so forth. So bp divided by 3 no remander should do the trick.
This will give the amino acid position or "aa_pos"
"""
data=pd.read_csv(snv)
df=data.sort_values("Id") # Sorting the data frame by sample name. NOTE : The indeces don't change
df["pos"]=df["coding.pos"]-1 # Note the pandas data frame is in base 1 while everthing else is in base 0 this fixes that
df["aa_pos"]=df["pos"]//3

"""
Here I am looping through the sample names (Id) in the variant csv. For each Id,
I am finding the index of that Id in the Id list; This is also the index of the sample in the protein list.
There are mulitple variants in some samples - by looping through the Ids in the data frame I will grab each of
these instances separately. I am using a counter to keep track of what row of the data frame. I am then taking
the appropriate amino acid position and grabing the consensus amino acid for that position in the sample.
These amino acids are stored in a list and then added to the data frame in the last line
"""

counter=0
consens=[]
for id in df.Id:
    if str(id) in Ids:
        i=Ids.index(str(id)) # get the index of the sample
        aa_pos=list(df["aa_pos"])[counter]
        aa=str(protein[i])[aa_pos]
        #print(str(i)+ " : " + str(id) + " : "+str(aa_pos) + " : " + aa)
        consens.append(aa)
    else:
        consens.append("NA")
    counter=counter+1
df["aa_consensus"]=consens


"""
Now I am looping through the variants by sample name (as above). This time I am introducing the mutation into the
nucleaic acid consensus. I am then retranslating the seqeunce, and grabing the appropriate amino acid. As before
each step is saved as a list and in the end the variant amino acid is added to the data frame.
If the variant amino acid is different than the consensus, the sample Id, nucleaic acid, and amino acid sequences
are also saved in separate lists to be used later (denoted as dn_*)
"""


counter=0
var_nucleotide=[]
var_Ids=[]
var_protein=[]
var_aa=[]

dn_nucleotide=[]
dn_Ids=[]
dn_protein=[]


for id in df.Id:
    if str(id) in Ids:
        i=Ids.index(str(id)) # get the index of the sample
        pos=list(df["pos"])[counter]
        var=list(df["var"])[counter]
        seq=list(str(nucleotide[i])) # get the consensus and make it a list so we can make the changes
        seq[pos]=var
        seq=''.join(seq)
        seq=Seq(seq,IUPAC.unambiguous_dna)
        var_nucleotide.append(seq)
        prot=seq.translate(to_stop=True)
        aa_pos=list(df["aa_pos"])
        aa_pos=aa_pos[counter]
        aa=str(prot)[aa_pos]
        con_aa=list(df["aa_consensus"])[counter]
        if aa !=con_aa:
            dn_protein.append(prot)
            dn_Ids.append(str(id))

        #print(str(id) + " : "+str(aa_pos) + " : " + aa + " : " + str(counter))
        var_aa.append(aa)
        var_protein.append(prot)
        var_Ids.append(str(id))
    else:
        var_nucleotide.append("NA")
        var_Ids.append("NA")
        var_protein.append("NA")
        var_aa.append("NA")
    counter=counter+1
df["aa_var"]=var_aa

"""
Id is the sample name, ref is the reference base in the plasmid control, var is the variant in the sample, freq.var
the frequency, pos the nucleaic acid postion, aa_pos the amino acid position, aa_consensus the amino acid at the
consensus level, and aa_var the variant amino acid.
"""

df[["Id","ref",'var',"freq.var","pos","aa_pos","aa_consensus","aa_var"]].to_csv("./results/"+args.prefix[0]+".aa.csv")

nonsense=df.loc[(df.aa_consensus!=df.aa_var)]





"""
Here I am adpapting the HA number script made by Jesse Bloom of the Bloom lab. Currently his script runs as an
exicutable that takes in a file containing variables - sequence, position of interest, and path to aligner.
I will adapt this into a function to be used here. I will also have to alter the output.
Currently I am just grabbing some of the text output as formating in his script and saving it in lists to be
added to the data frame.
"""



def Align(headers_seqs, progpath, program='PROBCONS', musclegapopen=None):
    """Performs a multiple sequence alignment of two or more sequences.

    By default, the protein sequences are aligned using PROBCONS.  This is
        probably the most accurate alignment program.  However, it is
        slow and consumes large amounts of memory if you are aligning
        a very large number of sequences (typically if you are aligning
        more than several hundred).  In that case, you may prefer to use
        MUSCLE instead.  You can choose between the two with the 'program'
        option.  If you decide to use MUSCLE, you can also align nucleotide
        sequences with this program.

    'headers_seqs' is a list specifying the names of the sequences that we
        want to align.  Each entry is a 2-tuple '(head, seq)' where 'head' is
        a header giving the sequence name and other information (might be empty)
        and 'seq' is a string giving the protein sequence.  The list must have
        at least 2 entries.

    'progpath' should specify a directory containing the alignment program executable,
        either PROBCONS or MUSCLE.  The PROBCONS executable is assumed to have
        the name "probcons" in this directory.  The MUSCLE executable is assumed to
        have the name "muscle" in this directory.

    'program' specifies what program to use for the alignment.  By default, it is
        "PROBCONS".  If you wish to use MUSCLE instead, set it to "MUSCLE".

    'musclegapopen' sets the MUSCLE gap openining penalty to the specified
        value. By default it is None, meaning we use the MUSCLE default penalty.
        You can also set it to a number; for example -100 will lead to fewer gaps.

    This executable is used to perform a multiple sequence alignment of the proteins
        with the default settings of either PROBCONS or MUSCLE.  The returned variable is a
        new list 'aligned_headers_seqs'.  Each entry is a 2-tuple '(head, aligned_seq)'.
        'head' has the same meaning as on input (the sequence header) and
        'aligned_seq' is the aligned sequence, with gaps inserted as '-'
        as appropriate.  Therefore, all of the 'aligned_seq' entries in
        'aligned_headers_seqs' are of the same length.  The entries in 'aligned_headers_seq'
        are in the same order as in the input list 'headers_seqs'.
    """
    if not (isinstance(headers_seqs, list) and len(headers_seqs) >= 2):
        raise ValueError('header_seqs does not specify a list with at least two entries.')
    if not os.path.isdir(progpath):
        raise ValueError("Cannot find directory %s." % progpath)
    if program == 'PROBCONS':
        exe = os.path.abspath("%s/probcons" % progpath) # the executable
    elif program == 'MUSCLE':
        exe = os.path.abspath("%s/muscle" % progpath) # the executable
    else:
        raise ValueError("Invalid value of %s for 'program'." % (str(program)))
    if not os.path.isfile(exe):
        raise IOError("Cannot find executable at %s." % exe)
    currdir = os.getcwd()
    tempdir = tempfile.mkdtemp()
    try:
        # do stuff in a temporary directory
        infile = "%s/in.fasta" % tempdir # input file
        outfile = "%s/out.fasta" % tempdir # output file
        WriteFASTA(headers_seqs, infile) # write sequences to the input file
        if program == 'PROBCONS':
            p = subprocess.Popen("%s %s" % (exe, infile), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # run ProbCons
            (output, errors) = p.communicate()

            open(outfile, 'w').write(output)
        elif program == 'MUSCLE':
            if musclegapopen != None:
                p = subprocess.Popen("%s -gapopen %d -in %s -out %s" % (exe, musclegapopen, infile, outfile), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # run MUSCLE
            else:
                p = subprocess.Popen("%s -in %s -out %s" % (exe, infile, outfile), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # run MUSCLE
            (output, errors) = p.communicate()
        try:
            aligned_headers_seqs = ReadFASTA(outfile)
        except:
            sys.stderr.write("Error getting alignment output, error of %s" % errors)
            raise
    finally:
        os.chdir(currdir) # return to the original directory
        for file in os.listdir(tempdir):
            os.remove("%s/%s" % (tempdir, file)) # remove files from temporary directory
        os.rmdir(tempdir) # remove temporary directory
    if len(aligned_headers_seqs) != len(headers_seqs):
        raise ValueError("Did not return the correct number of aligned sequences.")
    # put the aligned sequences in the same order as the input sequences
    n = len(aligned_headers_seqs[0][1]) # length of aligned sequences
    d = dict(aligned_headers_seqs)
    aligned_headers_seqs = []
    for (head, seq) in headers_seqs:
        try:
            alignedseq = d[head]
        except KeyError:
            raise ValueError("After alignment, the following header is missing: %s" % head)
        if len(alignedseq) != n:
            open('errors.temp', 'w').write(errors)
            raise ValueError("Aligned sequence %s is not of length %d: if you are using MUSCLE, you may be running out of memory.  Errors have been written to errors.temp." % (alignedseq, n))
        if len(seq) > n:
            open('errors.temp', 'w').write(errors)
            raise ValueError("Unaligned seq %s is longer than aligned length of %d: if you are using MUSCLE, you many be running out of memory.  Errors have been written to errors.temp." % (seq, n))
        aligned_headers_seqs.append((head, alignedseq))
    return aligned_headers_seqs # return the aligned sequences


def StripGapsToFirstSequence(aligned_headers_seqs):
    """Strips gaps from a reference sequence, and all corresponding alignments.

    On input, 'aligned_headers_seqs' should be a set of two or more aligned sequences,
        as would be returned by 'Align'.

    The first sequence in this alignment is taken to correspond to the reference sequence.
        The returned variable is a list similar to aligned_headers_seqs, but with
        all positions corresponding to gaps in this reference sequence stripped away.
        All gaps ('-') characters are removed from this reference sequence.  In addition,
        in all other aligned sequences in 'aligned_headers_seqs', every character at
        the same position as a gap in the reference sequence is removed.  Therefore,
        at the end of this procedure, all of the alignments have the same length
        as the reference sequence with its gaps stripped away.  The headers are
        unchanged.  The order of sequences in this stripped alignment is also
        unchanged.

    >>> StripGapsToFirstSequence([('s1', '-AT-A-GC'), ('s2', 'AAT-TAGC'), ('s3', '--T-A-GC')])
    [('s1', 'ATAGC'), ('s2', 'ATTGC'), ('s3', '-TAGC')]
    """
    if not (isinstance(aligned_headers_seqs, list) and len(aligned_headers_seqs) >= 2):
        raise ValueError("aligned_headers_seqs does not specify at least two aligned sequences.")
    (ref_head, ref_seq) = aligned_headers_seqs[0]
    non_strip_positions = [] # positions not to strip away
    stripped_ref_seq = []
    for i in range(len(ref_seq)):
        r = ref_seq[i]
        if r != '-':
            non_strip_positions.append(i)
            stripped_ref_seq.append(r)
    stripped_headers_seqs = [(ref_head, ''.join(stripped_ref_seq))]
    for (ihead, iseq) in aligned_headers_seqs[1 : ]:
        istrippedseq = ''.join([iseq[i] for i in non_strip_positions])
        stripped_headers_seqs.append((ihead, istrippedseq))
    return stripped_headers_seqs


def WriteFASTA(headers_seqs, filename, writable_file=False):
    """Writes sequences to a FASTA file.

    'headers_seqs' is a list of 2-tuples specifying sequences and their
        corresponding headers.  Each entry is the 2-tuple '(header, seq)'
        where 'header' is a string giving the header (without the leading ">"),
        and 'seq' is the corresponding sequence.

    'filename' is a string that specifies the name of the file to which the
        headers and sequences should be written.  If this file already exists,
        it is overwritten.

    'writable_file' is a Boolean switch specifying that rather than 'filename'
        giving a string specifying the name of a file to which the sequences
        should be written, it instead specifies a writable file object to which
        the sequences should be written.

    The sequences are written to the file in the same order that they are specified
        in 'headers_seqs'.
    """
    assert isinstance(writable_file, bool)
    if writable_file:
        f = filename
    else:
        f = open(filename, 'w')
    for (header, seq) in headers_seqs:
        f.write(">%s\n%s\n" % (header, seq))
    if not writable_file:
        f.close()



def ReadFASTA(fastafile):
    """Reads sequences from a FASTA file.

    'fastafile' should specify the name of a FASTA file.

    This function reads all sequences from the FASTA file.  It returns the
        list 'headers_seqs'.  This list is composed of a 2-tuple '(header, seq)'
        for every sequence entry in FASTA file.  'header' is the header for
        a sequence, with the leading ">" and any trailing spaces removed. 'seq'
        is the corresponding sequence.
    """
    lines = open(fastafile).readlines()
    headers_seqs = []
    header = None
    seq = []
    for line in lines:
        if line[0] == '>':
            if (not header) and (not seq):
                pass # first sequence in file
            elif header and not seq:
                raise ValueError("Empty sequence for %s" % header)
            elif seq and not header:
                raise ValueError("File does not begin with header.")
            else:
                seq = ''.join(seq)
                seq = seq.replace(' ', '')
                headers_seqs.append((header, seq))
            header = line.strip()[1 : ]
            seq = []
        else:
            seq.append(line.strip())
    if (not header) and (not seq):
        pass # first sequence in file
    elif header and not seq:
        raise ValueError("Empty sequence for %s" % header)
    elif seq and not header:
        raise ValueError("File does not begin with header.")
    else:
        seq = ''.join(seq)
        seq = seq.replace(' ', '')
        headers_seqs.append((header, seq))
    return headers_seqs


def GetCorrespondingResidue(seqs, i):
    """Gets the corresponding residue number for two aligned sequences.

    *seqs* is a set of two aligned sequences as *(head, sequence)* 2-tuples.

    *i* is the number of a residue in sequential numbering of *seqs[0]*
    without considering any of the gaps induced by alignment, in 1, 2, ...
    numbering.

    Returns the number of the residue in sequential numbering of *seqs[1]*
    without considering any of the gaps induced by alignment in 1, 2, ...
    numbering. Returns *None* if residue *i* of *seqs[0]* aligns with a
    gap in *seqs[1]*.
    """
    assert len(seqs) == 2
    s1 = seqs[0][1]
    s2 = seqs[1][1]
    assert len(s1) == len(s2)
    assert 1 <= i <= len(s1)
    s1index = s2index = 0
    for j in range(len(s1)):
        if s1[j] != '-':
            s1index += 1
        if s2[j] != '-':
            s2index += 1
        if s1index == i:
            if s2[j] == '-':
                return None
            else:
                return s2index


def find_aligments(aligner,aligner_path,seq,pos): # The positions will be given as base 0 and adjusted to match the convention (base 1) in the funciton
    """Main body of script."""
    print("\nBeginning execution of HA_numbering script.")

    # parse arguments
    if aligner == 'probcons':
        alignerpath = aligner_path # redundant
        if not os.path.isdir(alignerpath):
            raise IOError("The directory of %s specified by probconspath does not exist." % (alignerpath))
        prog = 'PROBCONS'
    elif aligner == 'muscle':
        alignerpath = aligner_path # redundant
        if not os.path.isdir(alignerpath):
            raise IOError("The directory of %s specified by musclepath does not exist." % (alignerpath))
        prog = 'MUSCLE'
    else:
        raise IOError("First argument does not specify probcons or muscle")
    ha_sequence = str(seq) # just get the sequence
    sites = pos
    # Define sequences and their numbering conversions.
    # The sequences are in seq_d and keyed by PDB code.
    # The number conversions from sequential numbering of these sequences to
    # th number labels are in label_d and keyed by PDB codes.
    seq_names = ['4HMG', '4JTV']
    seq_d = {}
    label_d = {}
    seq_4hmg_a = \
        'QDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQ' +\
        'NETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGST' +\
        'YPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGQSSRISIYWTIVKPG' +\
        'DVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGM' +\
        'RNVPEKQT'
    seq_4hmg_b = \
        'GLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDL' +\
        'EKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHD' +\
        'VYRDEALNNRFQIKG'
    seq_d['4HMG'] = seq_4hmg_a + seq_4hmg_b
    label_d['4HMG'] = dict([(i + 1, '%d in HA1' % (i + 1)) for i in range(len(seq_4hmg_a))] + [(len(seq_4hmg_a) + i + 1, '%d in HA2' % (i + 1)) for i in range(len(seq_4hmg_b))])
    assert len(seq_d['4HMG']) == len(label_d['4HMG'])
    seq_4jtv_a = \
        'DTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIV' +\
        'ETPSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPK' +\
        'LSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADTYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKI' +\
        'TFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRN' +\
        'I'
    seq_4jtv_b = \
        'GLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENL' +\
        'NKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYP' +\
        'KY'
    seq_d['4JTV'] = seq_4jtv_a + seq_4jtv_b
    label_d['4JTV'] = dict([(i + 1, '%d of HA1' % (i + 7)) for i in range(len(seq_4jtv_a))] + [(len(seq_4jtv_a) + i + 1, '%d of HA2' % (i + 1)) for i in range(len(seq_4jtv_b))])
    assert len(seq_d['4JTV']) == len(label_d['4JTV'])

    # make alignments
    print("Making %s alignments..." % prog)
    alignments = {}
    for seqname in seq_names:
        alignments[seqname] = Align([('seq', ha_sequence), (seqname, seq_d[seqname])], alignerpath, prog)
    print("Alignments complete.\n\nHere are the corresponding residue numbers:")
    for site in sites:
        site=site+1 # to correct for difference in base numbering
        if not (1 <= site <= len(ha_sequence)):
            raise ValueError("site %d is outside the valid range for sequential numbering of ha_sequence starting at 1." % site)
        sitestring = ['\nResidue %s%d in sequential numbering of your HA sequence corresponds to:' % (ha_sequence[site - 1], site)]
        for seqname in seq_names:
            i = GetCorrespondingResidue(alignments[seqname], site)
            if i == None:
                sitestring.append('  * an alignment gap in %s' % seqname)
            else:
                if not (1 <= i <= len(seq_d[seqname])):
                    raise ValueError("Invalid corrresponding residue for %s -- something is wrong with this program" % seqname)
                aa = seq_d[seqname][i - 1]

                sitestring.append('  * %s%s in %s' % (aa, label_d[seqname][i], seqname))
        print( '\n'.join(sitestring))
        return([sitestring[1][4:-7],sitestring[2][4:-7]])

"""
Now I am looping through all nonsynonymous mutations as before and applying the function above to get the
H3 numbering system. The script takes in base 0 positions and handels them accordingly.
"""


counter=0
PDB_4HMG=[]
PDB_4JTV=[]
for id in nonsense.Id:
    if str(id) in dn_Ids:
        i=dn_Ids.index(str(id)) # get the index of the sample
        pos=list(nonsense["aa_pos"])[counter]
        nice=find_aligments('muscle',muscle_dir,dn_protein[i],[pos])
        PDB_4HMG.append(nice[0])
        PDB_4JTV.append(nice[1])
    counter=counter+1
nonsense["PDB_4HMG"] = PDB_4HMG
nonsense["PDB_4JTV"]=PDB_4JTV

print(nonsense[["Id","ref",'var',"freq.var","pos","aa_pos","aa_consensus","aa_var","PDB_4HMG","PDB_4JTV"]])

nonsense[["Id","ref",'var',"freq.var","pos","aa_pos","aa_consensus","aa_var","PDB_4HMG","PDB_4JTV"]].to_csv("./results/"+args.prefix[0]+".putative.antigenic.csv")
