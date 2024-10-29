#!/usr/bin/sh
INPUT_DIR="/home/shuaib/L4T1/hg002_benchmark_files"
OUTPUT_DIR="${PWD}/Outputs"

TRUTH_BED=${INPUT_DIR}/chr1_22_full.bed
REF=${INPUT_DIR}/hs37d5.fa

# OUTPUT DIRECTORIES WHERE TRAINING EXAMPLES WILL BE SAVED
TRAIN_OUTPUT=${INPUT_DIR}/PEPPER_TRAIN_IMAGES/final_runs/tmp
# TEST_OUTPUT=PEPPER_TEST_IMAGES/cf_0175
TEST_OUTPUT=PEPPER_TEST_IMAGES/

THREADS=1

for coverage in 50
do
  BAM=${INPUT_DIR}/HG002_35x_HiFi_2_GRCh37.bam

  echo ${BAM}

  time pepper_variant make_images \
  -b ${BAM} \
  -f ${REF} \
  -r 15:54201236-54205127 \
  -rb ${TRUTH_BED} \
  -t ${THREADS} \
  -o ${TRAIN_OUTPUT} \
  -d 1.0 \
  -p 1.0 \
  --hifi
done

  # --include_supplementary