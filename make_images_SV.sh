#!/usr/bin/sh
INPUT_DIR="/home/shuaib/L4T1/hg002_benchmark_files"
OUTPUT_DIR="${PWD}/Outputs"

TRUTH_BED=${INPUT_DIR}/SVs_only_passed.bed
TRUTH_VCF=${INPUT_DIR}/SVs_only_passed.vcf.gz
REF=${INPUT_DIR}/hs37d5.fa

# OUTPUT DIRECTORIES WHERE TRAINING EXAMPLES WILL BE SAVED
TRAIN_OUTPUT=${INPUT_DIR}/PEPPER_TRAIN_IMAGES/high_conf_only
TEST_OUTPUT=${INPUT_DIR}/PEPPER_TEST_IMAGES/termbreak

THREADS=1

for coverage in 50
do
  BAM=${INPUT_DIR}/HG002_35x_HiFi_2_GRCh37.bam

  echo ${BAM}

  time pepper_variant_train make_train_images \
  -b ${BAM} \
  -f ${REF} \
  -tv ${TRUTH_VCF} \
  -r 4:108276163-108280199 \
  -rb ${TRUTH_BED} \
  -t ${THREADS} \
  -o ${TRAIN_OUTPUT} \
  -d 1.0 \
  -p 1.0 \
  --hifi  
done
