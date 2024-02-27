#!/usr/bin/sh
INPUT_DIR="/home/fahmid/Thesis/Data"
OUTPUT_DIR="/home/fahmid/Thesis/Output"

BAM_DIR="/mnt/f"

# TRUTH_BED=${INPUT_DIR}/SVs_only_passed.bed
TRUTH_BED=${INPUT_DIR}/HG002_SVs_Tier1_v0.6.bed
TRUTH_VCF=${INPUT_DIR}/SVs_only_passed.vcf.gz
REF=${INPUT_DIR}/hs37d5.fa

# OUTPUT DIRECTORIES WHERE TRAINING EXAMPLES WILL BE SAVED
TRAIN_OUTPUT=${OUTPUT_DIR}/PEPPER_TRAIN_IMAGES
TEST_OUTPUT=${OUTPUT_DIR}/PEPPER_TEST_IMAGES

THREADS=16

for coverage in 50
do
  BAM=${BAM_DIR}/hg002_74gb/HG002_35x_HiFi_2_GRCh37.bam

  echo ${BAM}

  time pepper_variant make_images \
  -b ${BAM} \
  -f ${REF} \
  -r 1:75837899-75851590\
  -rb ${TRUTH_BED} \
  -t ${THREADS} \
  -o ${TRAIN_OUTPUT} \
  -d 1.0 \
  -p 1.0 \
  --hifi  
done
