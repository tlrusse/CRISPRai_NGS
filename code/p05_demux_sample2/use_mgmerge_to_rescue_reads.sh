set -euo pipefail

mkdir -p ${DATA_OUT}/15754-PP/
FASTQ1=${DATA_IN}/fastqs_15754-PP/15754-PP-1_S1_R1_001.fastq.gz
FASTQ2=${DATA_IN}/fastqs_15754-PP/15754-PP-1_S1_R2_001.fastq.gz
OUT_FASTQ=${DATA_OUT}/15754-PP/15754-PP-1_S1_001_pe_ngmerged.fastq.gz

# because this run has some reads with low quality, we might be able to 
# rescue some of those reads by merging the paired-end reads together, 
# statistically based on phred scores.

# note: for this run, we have 150bp paired-end reads that are trying to capture a 280 bp amplicon, so we expect a 20bp overlap between the reads.
NGmerge -1 ${FASTQ1} -2 ${FASTQ2} -o ${OUT_FASTQ} -d -e 20 -n 64 -v

whecho "finished merging reads with NGmerge, output is ${OUT_FASTQ}"