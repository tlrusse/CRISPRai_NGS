set -euo pipefail

mkdir -p ${DATA_OUT}/15107-PP/
mkdir -p ${DATA_OUT}/15107-PP/alignment/
FASTQ1=${DATA_IN}/fastqs_15107-PP/15107-PP-1_S1_R1_001.fastq.gz
FASTQ2=${DATA_IN}/fastqs_15107-PP/15107-PP-1_S1_R2_001.fastq.gz
OUT_SAM=${DATA_OUT}/15107-PP/alignment/15107-PP-1_S1_001_pe.sam
OUT_LOG=${DATA_OUT}/15107-PP/alignment/15107-PP-1_S1_001_pe_bowtie2.log

# use for unpaired reads
# bowtie2 -p 32 -x barcode_primer_ind -U ${DATA_OUT}9536-PP_reseq_processing/flow_cell_chunks/9536-PP-1__S1_R1_001_flow_cell_${FLOW_CELL}.fastq --local --score-min L,30,0 > ${OUT_SAM} 2> ${OUT_LOG}
# use for paired-end reads (note that the sam file will have mappings for both reads)
bowtie2 -p 32 -x barcode_primer_ind -1 ${FASTQ1} -2 ${FASTQ2} --local --score-min L,30,0 > ${OUT_SAM} 2> ${OUT_LOG}

samtools stats ${OUT_SAM} | grep ^SN > ${OUT_SAM}.summary