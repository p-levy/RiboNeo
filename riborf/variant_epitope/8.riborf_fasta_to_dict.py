#!/usr/bin/env python3

from Bio import SeqIO
import json

# AA dict
aa_dict = dict()
with open("/Volumes/labagros/pierre/proj/riboseq/riborf/360RIO136_gencode_out/translating_candidateORF_AA.fa") as input_file:
	my_list = list(SeqIO.parse(input_file, "fasta"))

	for record in my_list:
		new_record = record.id.replace(":", "|")
		aa_dict[new_record] = str(record.seq)

	json.dump(aa_dict, open("/Volumes/labagros/pierre/proj/riboseq/riborf/360RIO136_gencode_out/translating_candidateORF_AA.dict",'w'))

# CDS dict
cds_dict = dict()
with open("/Volumes/labagros/pierre/proj/riboseq/riborf/360RIO136_gencode_out/translating_candidateORF.fa") as input_file:
	my_list = list(SeqIO.parse(input_file, "fasta"))

	for record in my_list:
		new_record = record.id.replace(":", "|")
		cds_dict[new_record] = str(record.seq)

	json.dump(cds_dict, open("/Volumes/labagros/pierre/proj/riboseq/riborf/360RIO136_gencode_out/translating_candidateORF.dict",'w'))
