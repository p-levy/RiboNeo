# Build Gedi/PRICE Docker Image
`cd docker`

`docker build -t docker_image:tag .` (replace docker_image:tag as desired)

# Prepare Genome

Done only once per genome. Replace paths and docker image name in following command.

```
docker run --rm \
	-v /local/path/to/refs:/refs \
	docker_image:tag \
	gedi -e IndexGenome -s /refs/path/to/genome.fasta \
    -a /refs/path/to/annotation.gtf \
    -n prefixOutput \
    -o /refs/prefixOutput.oml \
    -nobowtie \
    -nostar \
    -nokallisto \
    -p
```

# Run PRICE

 Replace paths and docker image name in following command.

``` 
docker run --rm \
    -v /path/to/bam/dir:/bamdir \
	-v /path/to/output/dir:/outdir \
	-v /local/path/to/refs:/refs \
	-w /outdir \
	docker_image:tag \
	gedi -e Price -reads /bamdir/Aligned.sortedByCoord.out.bam \
	-genomic /refs/prefixOutput.oml \
	-prefix outputPrefix \
	-progress -plot
```

# PRICE post-process
## Select ORFs with FDR <10%
Use ViewCIT utility of gedi to extract list of orfs passing the FDR cutoff of 10% (default) from the orfs.cit file.
```
OUTDIR=/path/to/output/
docker run --rm \
	-v $OUTDIR:/outdir \
	-w /outdir \
	docker_image:tag \
    /bin/bash -c "gedi -e ViewCIT -m bed -name 'd.getType()'  outputPrefix.orfs.cit > outputPrefix.orfs.fdr10.bed"
```

## Select noncanonical ORFs (NOTE: not needed maybe)
egrep -v 'CDS|Ext|Trunc|Variant|iORF' outputPrefix.orfs.fdr10.bed > outputPrefix.NCorfs.fdr10.bed

