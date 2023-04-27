# If seqkit is to be ran with docker

docker run --rm -v /mnt/bioinfnas/immuno/plevy/proj/riboseq/riborf/360RIO136_gencode_out:/wdir \
	pegi3s/seqkit:latest grep -f /wdir/translating_ORF_ID.txt /wdir/candidateORF.fa \
	> /mnt/bioinfnas/immuno/plevy/proj/riboseq/riborf/360RIO136_gencode_out/translating_candidateORF.fa