from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


for seq_record in SeqIO.parse("/Users/jt/Documents/Analysis/FluVac_paper/data/processed/Run_1293/parsed_fa/HA.fa", "fasta"):
    print(seq_record.id)
    seq_record.seq.alphabet=IUPAC.unambiguous_dna
    fixed = list(str(seq_record.seq)) # the sequences all start as ATT they should be ATG - this is likely due to issues of coverage near the end. Fixing this error opens up the proper reading frame and when done to the plasmid control it yeilds a 100% match across all 566aa found in NCBI
    fixed[2]="G"
    fixed=''.join(fixed)
    fixed_record=Seq(fixed,IUPAC.unambiguous_dna)
    prot=fixed_record.translate()
    print prot
    print(len(prot))
