# make sure to install mageck into your conda environment first
# conda install -c bioconda -c conda-forge mageck

LIBRARY_FILE=${DATA_OUT}/15107-PP/library_files/ # finish this file path

BC1_READS=${DATA_OUT}/15107-PP/alignment/unpaired_demux_R1_only/15107_PP_sample_demux_barcode_ACAGTG.fastq
BC2_READS=${DATA_OUT}/15107-PP/alignment/unpaired_demux_R1_only/15107_PP_sample_demux_barcode_CGATGT.fastq
BC3_READS=${DATA_OUT}/15107-PP/alignment/unpaired_demux_R1_only/15107_PP_sample_demux_barcode_GCCAAT.fastq
BC4_READS=${DATA_OUT}/15107-PP/alignment/unpaired_demux_R1_only/15107_PP_sample_demux_barcode_TGACCA.fastq

BC1_NAME="bc1"
BC2_NAME="bc2"
BC3_NAME="bc3"
BC4_NAME="bc4"

mkdir -p ${DATA_OUT}/mageck/

# if barcode 1 and 2 came from the same sample, use the following mageck commands
# mageck count -l ${LIBRARY_FILE}/15107-PP_library_file.txt -n ${DATA_OUT}/mageck/15107_PP_${BC1_NAME}_${BC2_NAME} --fastq ${BC1_READS},${BC2_READS} --sample-label ${BC1_NAME},${BC2_NAME}

# if barcodes 1 and 2 came from different samples, use the following mageck commands
# mageck count -l ${LIBRARY_FILE}/15107-PP_library_file.txt -n ${DATA_OUT}/mageck/15107_PP_${BC1_NAME} --fastq ${BC1_READS} --sample-label ${BC1_NAME}