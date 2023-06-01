import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Create amino acid fasta from nucleotide fasta for translating ORFs',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('outdir',
                        metavar='outdir',
                        help='Path to output directory')

    return parser.parse_args()

# --------------------------------------------------

# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    outdir = args.outdir

    with open(f"{outdir}/translating_candidateORF_AA.fa", 'w') as aa_fa:
        for dna_record in SeqIO.parse(f"{outdir}/translating_candidateORF.fa", 'fasta'):
            dna_seq = dna_record.seq

            # generate translation
            aa_seq = dna_seq.translate()

            # write new record
            aa_record = SeqRecord(aa_seq, id=dna_record.id, description="")
            SeqIO.write(aa_record, aa_fa, 'fasta')

# --------------------------------------------------
if __name__ == '__main__':
    main()
