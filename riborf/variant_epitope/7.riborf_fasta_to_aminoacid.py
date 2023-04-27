import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

with open("/Volumes/labagros/pierre/proj/riboseq/riborf/360RIO136_gencode_out/translating_candidateORF_AA.fa", 'w') as aa_fa:
    for dna_record in SeqIO.parse("/Volumes/labagros/pierre/proj/riboseq/riborf/360RIO136_gencode_out/translating_candidateORF.fa", 'fasta'):
        dna_seq = dna_record.seq

        # generate translation
        aa_seq = dna_seq.translate()

        # write new record
        aa_record = SeqRecord(aa_seq, id=dna_record.id, description="")
        SeqIO.write(aa_record, aa_fa, 'fasta')
