#!/usr/bin/sh
INPUT_DIR="/home/fahmid/Thesis/Data"
OUTPUT_DIR="/home/fahmid/Thesis/Output"

BAM_DIR="/mnt/f"

TRUTH_BED=${INPUT_DIR}/SVs_only_passed.bed
TRUTH_VCF=${INPUT_DIR}/SVs_only_passed.vcf.gz
REF=${INPUT_DIR}/hs37d5.fa

# OUTPUT DIRECTORIES WHERE TRAINING EXAMPLES WILL BE SAVED
TRAIN_OUTPUT=${OUTPUT_DIR}/PEPPER_TRAIN_IMAGES
TEST_OUTPUT=${OUTPUT_DIR}/PEPPER_TEST_IMAGES

THREADS=1

for coverage in 50
do
  BAM=${BAM_DIR}/hg002_74gb/HG002_35x_HiFi_2_GRCh37.bam

  echo ${TRUTH_BED}
  echo ${TRUTH_VCF}
  echo ${REF}
  echo ${BAM}

  time pepper_variant_train make_train_images \
  -b ${BAM} \
  -f ${REF} \
  -tv ${TRUTH_VCF} \
  -rb ${TRUTH_BED} \
  -r 1:86740185-86742931 \
  -t ${THREADS} \
  -o ${TEST_OUTPUT} \
  -d 1.0 \
  -p 1.0 \
  --hifi  
done

  # -r 1:58743094-58745652\
