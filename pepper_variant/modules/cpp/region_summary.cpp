//
// Created by Kishwar Shafin on 3/21/21.
//

#include "region_summary.h"

#include <utility>

RegionalSummaryGenerator::RegionalSummaryGenerator(long long region_start, long long region_end, string reference_sequence) {
    this->ref_start = region_start;
    this->ref_end = region_end;
    this->reference_sequence = std::move(reference_sequence);
    this->total_observered_insert_bases = 0;
    this->max_observed_insert.resize(region_end-region_start+1, 0);
    this->cumulative_observed_insert.resize(region_end-region_start+1, 0);
}

void RegionalSummaryGenerator::generate_max_insert_observed(const type_read& read) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    long long reference_index;
    if(ImageOptionsRegion::GENERATE_INDELS == true) {
        for (auto &cigar: read.cigar_tuples) {
            if (ref_position > ref_end) break;
            switch (cigar.operation) {
                case CIGAR_OPERATIONS::EQUAL:
                case CIGAR_OPERATIONS::DIFF:
                case CIGAR_OPERATIONS::MATCH:
                    cigar_index = 0;
                    if (ref_position < ref_start) {
                        cigar_index = min(ref_start - ref_position, (long long) cigar.length);
                        read_index += cigar_index;
                        ref_position += cigar_index;
                    }
                    for (int i = cigar_index; i < cigar.length; i++) {
                        reference_index = ref_position - ref_start;
                        read_index += 1;
                        ref_position += 1;
                    }
                    break;
                case CIGAR_OPERATIONS::IN:
                    reference_index = ref_position - ref_start - 1;

                    if (ref_position - 1 >= ref_start &&
                        ref_position - 1 <= ref_end) {
                        max_observed_insert[reference_index] = std::max(max_observed_insert[reference_index], (uint64_t) cigar.length);
                    }
                    read_index += cigar.length;
                    break;
                case CIGAR_OPERATIONS::REF_SKIP:
                case CIGAR_OPERATIONS::PAD:
                case CIGAR_OPERATIONS::DEL:
                    reference_index = ref_position - ref_start - 1;
                    ref_position += cigar.length;
                    break;
                case CIGAR_OPERATIONS::SOFT_CLIP:
                    read_index += cigar.length;
                    break;
                case CIGAR_OPERATIONS::HARD_CLIP:
                    break;
            }
        }
    }
}


void RegionalSummaryGenerator::generate_max_insert_summary(vector <type_read> &reads) {
    for (auto &read:reads) {
        // this populates base_summaries and insert_summaries dictionaries
        generate_max_insert_observed(read);
    }

    cumulative_observed_insert[0] = 0;
    total_observered_insert_bases += max_observed_insert[0];

    positions.push_back(ref_start);
    index.push_back(0);

    for(int j=1; j <= max_observed_insert[0]; j++) {
        positions.push_back(ref_start);
        index.push_back(j);
    }

    for(int i=1;i < max_observed_insert.size(); i++) {
        cumulative_observed_insert[i] = cumulative_observed_insert[i-1] + max_observed_insert[i-1];
        total_observered_insert_bases += max_observed_insert[i];
        positions.push_back(ref_start + i);
        index.push_back(0);
        for(int j=1; j <= max_observed_insert[i]; j++) {
            positions.push_back(ref_start + i);
            index.push_back(j);
        }
    }
}


char check_truth_base(char base) {
    if(base=='A' || base=='a' ||
       base=='C' || base=='c' ||
       base=='T' || base=='t' ||
       base=='G' || base=='g' ||
       base=='*' || base=='#') return base;
    return '*';
}

uint8_t get_label_index(char base_h1, char base_h2) {
    base_h1 = toupper(base_h1);
    base_h2 = toupper(base_h2);
    // AA
    if (base_h1 == 'A' && base_h2 == 'A') return 1;
    // AC
    if (base_h1 == 'A' && base_h2 == 'C') return 2;
    if (base_h1 == 'C' && base_h2 == 'A') return 2;
    // AT
    if (base_h1 == 'A' && base_h2 == 'T') return 3;
    if (base_h1 == 'T' && base_h2 == 'A') return 3;
    // AG
    if (base_h1 == 'A' && base_h2 == 'G') return 4;
    if (base_h1 == 'G' && base_h2 == 'A') return 4;
    // A*
    if (base_h1 == 'A' && base_h2 == '*') return 5;
    if (base_h1 == '*' && base_h2 == 'A') return 5;
    // A#
    if (base_h1 == 'A' && base_h2 == '#') return 6;
    if (base_h1 == '#' && base_h2 == 'A') return 6;
    // CC
    if (base_h1 == 'C' && base_h2 == 'C') return 7;
    // CT
    if (base_h1 == 'C' && base_h2 == 'T') return 8;
    if (base_h1 == 'T' && base_h2 == 'C') return 8;
    // CG
    if (base_h1 == 'C' && base_h2 == 'G') return 9;
    if (base_h1 == 'G' && base_h2 == 'C') return 9;
    // C*
    if (base_h1 == 'C' && base_h2 == '*') return 10;
    if (base_h1 == '*' && base_h2 == 'C') return 10;
    // C#
    if (base_h1 == 'C' && base_h2 == '#') return 11;
    if (base_h1 == '#' && base_h2 == 'C') return 11;
    // TT
    if (base_h1 == 'T' && base_h2 == 'T') return 12;
    // TG
    if (base_h1 == 'T' && base_h2 == 'G') return 13;
    if (base_h1 == 'G' && base_h2 == 'T') return 13;
    // T*
    if (base_h1 == 'T' && base_h2 == '*') return 14;
    if (base_h1 == '*' && base_h2 == 'T') return 14;
    // T#
    if (base_h1 == 'T' && base_h2 == '#') return 15;
    if (base_h1 == '#' && base_h2 == 'T') return 15;
    // GG
    if (base_h1 == 'G' && base_h2 == 'G') return 16;
    // G*
    if (base_h1 == 'G' && base_h2 == '*') return 17;
    if (base_h1 == '*' && base_h2 == 'G') return 17;
    // G#
    if (base_h1 == 'G' && base_h2 == '#') return 18;
    if (base_h1 == '#' && base_h2 == 'G') return 18;
    // **
    if (base_h1 == '*' && base_h2 == '*') return 19;
    // *#
    if (base_h1 == '*' && base_h2 == '#') return 20;
    if (base_h1 == '#' && base_h2 == '*') return 20;
    // ##
    if (base_h1 == '#' && base_h2 == '#') return 0;
    return 0;
}


uint8_t get_variant_type_label_index(int type_h1, int type_h2) {

    if (type_h1 == VariantTypes::HOM_REF && type_h2 == VariantTypes::HOM_REF) return 0;

    if (type_h1 == VariantTypes::HOM_REF && type_h2 == VariantTypes::SNP) return 1;
    if (type_h1 == VariantTypes::SNP && type_h2 == VariantTypes::HOM_REF) return 1;

    if (type_h1 == VariantTypes::HOM_REF && type_h2 == VariantTypes::INSERT) return 2;
    if (type_h1 == VariantTypes::INSERT && type_h2 == VariantTypes::HOM_REF) return 2;

    if (type_h1 == VariantTypes::HOM_REF && type_h2 == VariantTypes::DELETE) return 3;
    if (type_h1 == VariantTypes::DELETE && type_h2 == VariantTypes::HOM_REF) return 3;

    if (type_h1 == VariantTypes::SNP && type_h2 == VariantTypes::SNP) return 4;

    if (type_h1 == VariantTypes::SNP && type_h2 == VariantTypes::INSERT) return 5;
    if (type_h1 == VariantTypes::INSERT && type_h2 == VariantTypes::SNP) return 5;

    if (type_h1 == VariantTypes::SNP && type_h2 == VariantTypes::DELETE) return 6;
    if (type_h1 == VariantTypes::DELETE && type_h2 == VariantTypes::SNP) return 6;

    if (type_h1 == VariantTypes::INSERT && type_h2 == VariantTypes::INSERT) return 7;

    if (type_h1 == VariantTypes::INSERT && type_h2 == VariantTypes::DELETE) return 8;
    if (type_h1 == VariantTypes::DELETE && type_h2 == VariantTypes::INSERT) return 8;

    if (type_h1 == VariantTypes::DELETE && type_h2 == VariantTypes::DELETE) return 9;

    cout<<"ERROR: VARIANT LABEL NOT DEFINED: "<<type_h1<<" "<<type_h2<<endl;
    exit(1);
}


int RegionalSummaryGenerator::get_reference_feature_index(char base) {
    base = toupper(base);
    if (base == 'A') return ImageOptionsRegion::REFERENCE_INDEX_START;
    if (base == 'C') return ImageOptionsRegion::REFERENCE_INDEX_START + 1;
    if (base == 'G') return ImageOptionsRegion::REFERENCE_INDEX_START + 2;
    if (base == 'T') return ImageOptionsRegion::REFERENCE_INDEX_START + 3;
    return ImageOptionsRegion::REFERENCE_INDEX_START + 4;
}

void RegionalSummaryGenerator::encode_reference_bases(int **image_matrix) {
    for (long long ref_position = ref_start; ref_position <= ref_end; ref_position++) {
        // encode the C base
        int base_index = (int) (ref_position - ref_start + cumulative_observed_insert[ref_position - ref_start]);
        int feature_index = get_reference_feature_index(reference_sequence[ref_position - ref_start]);
        image_matrix[base_index][feature_index] = 255;

        for(int i = 1; i <= max_observed_insert[ref_position - ref_start]; i++) {
            base_index = (int) (ref_position - ref_start + cumulative_observed_insert[ref_position - ref_start]) + i;
            feature_index = get_reference_feature_index('*');
            image_matrix[base_index][feature_index] = 255;
        }
    }
}

int RegionalSummaryGenerator::get_feature_index(char base, bool is_reverse) {
    base = toupper(base);
    if (is_reverse) {
        if (base == 'A') return ImageOptionsRegion::BASE_INDEX_START;
        if (base == 'C') return ImageOptionsRegion::BASE_INDEX_START + 1;
        if (base == 'G') return ImageOptionsRegion::BASE_INDEX_START + 2;
        if (base == 'T') return ImageOptionsRegion::BASE_INDEX_START + 3;
        if (base == 'I') return ImageOptionsRegion::BASE_INDEX_START + 5;
        if (base == 'D') return ImageOptionsRegion::BASE_INDEX_START + 6;
        return ImageOptionsRegion::BASE_INDEX_START + 4;
    } else {
        // tagged and forward
        if (base == 'A') return ImageOptionsRegion::BASE_INDEX_START + 7;
        if (base == 'C') return ImageOptionsRegion::BASE_INDEX_START + 8;
        if (base == 'G') return ImageOptionsRegion::BASE_INDEX_START + 9;
        if (base == 'T') return ImageOptionsRegion::BASE_INDEX_START + 10;
        if (base == 'I') return ImageOptionsRegion::BASE_INDEX_START + 12;
        if (base == 'D') return ImageOptionsRegion::BASE_INDEX_START + 13;
        return ImageOptionsRegion::BASE_INDEX_START + 11;
    }
}



void RegionalSummaryGenerator::generate_labels(const vector<type_truth_record>& hap1_records, const vector<type_truth_record>& hap2_records) {
    int region_size = (int) (ref_end - ref_start + total_observered_insert_bases + 1);
    labels_hp1.resize(region_size + 1, '*');
    labels_hp2.resize(region_size + 1, '*');
    variant_type_labels_hp1.resize(region_size + 1, VariantTypes::HOM_REF);
    variant_type_labels_hp2.resize(region_size + 1, VariantTypes::HOM_REF);

    for(long long pos = ref_start; pos <= ref_end; pos++) {
        int base_index = (int)(pos - ref_start + cumulative_observed_insert[pos - ref_start]);
        labels_hp1[base_index] = reference_sequence[base_index];
        labels_hp2[base_index] = reference_sequence[base_index];
    }

    for(const auto& truth_record : hap1_records) {
        if(truth_record.ref.length() > truth_record.alt.length()) {
            //it's a delete
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp1[base_index] = VariantTypes::DELETE;
                labels_hp1[base_index] = '#';
            }
            // dont expand to the rest of the bases

//            for(long long pos = truth_record.pos_start; pos < truth_record.pos_end; pos++) {
//                if (pos >= ref_start && pos <= ref_end) {
//                    int base_index = (int) (pos - ref_start + cumulative_observed_insert[pos - ref_start]);
//                    if (pos - truth_record.pos_start < truth_record.alt.length()) {
//                        labels_hp1[base_index] = truth_record.alt[pos - truth_record.pos_start];
//                    } else {
//                        labels_hp1[base_index] = '*';
//                    }
//                }
//            }
        } else if(truth_record.ref.length() < truth_record.alt.length()) {
            //it's an insert
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp1[base_index] = VariantTypes::INSERT;
                labels_hp1[base_index] = '*';
            }

//            for(long long pos = truth_record.pos_start; pos < truth_record.pos_end; pos++) {
//                if (pos >= ref_start && pos <= ref_end) {
//                    int base_index = (int) (pos - ref_start + cumulative_observed_insert[pos - ref_start]);
//                    labels_hp1[base_index] = truth_record.alt[pos - truth_record.pos_start];
//                }
//            }
        } else if(truth_record.ref.length() == truth_record.alt.length()) {
            //it's a SNP
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp1[base_index] = VariantTypes::SNP;
            }

            for(long long pos = truth_record.pos_start; pos < truth_record.pos_end; pos++) {
                if (pos >= ref_start && pos <= ref_end) {
                    int base_index = (int) (pos - ref_start + cumulative_observed_insert[pos - ref_start]);
                    labels_hp1[base_index] = truth_record.alt[pos - truth_record.pos_start];
                }
            }
        }
    }

    for(const auto& truth_record : hap2_records) {
        if(truth_record.ref.length() > truth_record.alt.length()) {
            //it's a delete
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp2[base_index] = VariantTypes::DELETE;
                labels_hp2[base_index] = '#';
            }
            // dont expand to the rest of the bases

//            for(long long pos = truth_record.pos_start; pos < truth_record.pos_end; pos++) {
//                if (pos >= ref_start && pos <= ref_end) {
//                    int base_index = (int) (pos - ref_start + cumulative_observed_insert[pos - ref_start]);
//                    if (pos - truth_record.pos_start < truth_record.alt.length()) {
//                        labels_hp2[base_index] = truth_record.alt[pos - truth_record.pos_start];
//                    } else {
//                        labels_hp2[base_index] = '*';
//                    }
//                }
//            }
        } else if(truth_record.ref.length() < truth_record.alt.length()) {
            //it's an insert
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp2[base_index] = VariantTypes::INSERT;
                labels_hp2[base_index] = '*';
            }

//            for(long long pos = truth_record.pos_start; pos < truth_record.pos_end; pos++) {
//                if (pos >= ref_start && pos <= ref_end) {
//                    int base_index = (int) (pos - ref_start + cumulative_observed_insert[pos - ref_start]);
//                    labels_hp2[base_index] = truth_record.alt[pos - truth_record.pos_start];
//                }
//            }
        } else if(truth_record.ref.length() == truth_record.alt.length()) {
            //it's a SNP
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp2[base_index] = VariantTypes::SNP;
            }

            for(long long pos = truth_record.pos_start; pos < truth_record.pos_end; pos++) {
                if (pos >= ref_start && pos <= ref_end) {
                    int base_index = (int) (pos - ref_start + cumulative_observed_insert[pos - ref_start]);
                    labels_hp2[base_index] = truth_record.alt[pos - truth_record.pos_start];
                }
            }
        }
    }
}


void RegionalSummaryGenerator::populate_summary_matrix(int **image_matrix,
                                                       int *coverage_vector,
                                                       type_read read) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    double base_quality = 1.0 - pow(10.0 , (-1.0 * read.base_qualities[read_index])/10);
    for (int cigar_i=0; cigar_i<read.cigar_tuples.size(); cigar_i++) {
        CigarOp cigar = read.cigar_tuples[cigar_i];
        if (ref_position > ref_end) break;
        switch (cigar.operation) {
            case CIGAR_OPERATIONS::EQUAL:
            case CIGAR_OPERATIONS::DIFF:
            case CIGAR_OPERATIONS::MATCH:
                cigar_index = 0;
                if (ref_position < ref_start) {
                    cigar_index = min(ref_start - ref_position, (long long) cigar.length);
                    read_index += cigar_index;
                    ref_position += cigar_index;
                }
                for (int i = cigar_index; i < cigar.length; i++) {
                    //read.base_qualities[read_index] base quality
//                    double base_quality = 1.0 - pow(10.0 , (-1.0 * read.base_qualities[read_index])/10);

                    if (ref_position >= ref_start && ref_position <= ref_end) {
                        char base = read.sequence[read_index];
                        char ref_base = reference_sequence[ref_position - ref_start];

                        int base_index = (int)(ref_position - ref_start + cumulative_observed_insert[ref_position - ref_start]);
                        int feature_index = get_feature_index(base, read.flags.is_reverse);

                        //look forward and make sure this is not an anchor base
                        bool check_this_base = true;
                        if(i == cigar.length - 1 && cigar_i + 1 < read.cigar_tuples.size()) {
                            CigarOp next_cigar = read.cigar_tuples[cigar_i + 1];
                            if(next_cigar.operation == CIGAR_OPERATIONS::IN ||
                               next_cigar.operation == CIGAR_OPERATIONS::DEL) {
                                // this is an anchor base of a delete or an insert, don't process this.
                                coverage_vector[ref_position - ref_start] += 1;
                                check_this_base = false;
                            }
                        }

                        if(check_this_base) {
                            // update the summary of base
                            if (ref_position >= ref_start && ref_position <= ref_end) {
                                image_matrix[base_index][feature_index] += 1;
                                coverage_vector[ref_position - ref_start] += 1;
                            }
                        }

                    }
                    read_index += 1;
                    ref_position += 1;
                }
                break;
            case CIGAR_OPERATIONS::IN:
                if (ref_position - 1 >= ref_start && ref_position - 1 <= ref_end) {
                    // process insert allele here
                    string alt;
                    char ref_base = reference_sequence[ref_position - 1 - ref_start];
                    int base_index = (int)((ref_position - 1) - ref_start + cumulative_observed_insert[(ref_position - 1) - ref_start]);
                    int insert_count_index =  get_feature_index('I', read.flags.is_reverse);
                    image_matrix[base_index][insert_count_index] += 1;

//                    if(ImageOptionsRegion::GENERATE_INDELS == true) {
//                        alt = read.sequence.substr(read_index, cigar.length);
//                        for (int i = 0; i < cigar.length; i++) {
//                            base_index = (int) ((ref_position - 1) - ref_start + cumulative_observed_insert[(ref_position - 1) - ref_start] + (i + 1));
//                            int feature_index = get_feature_index(alt[i], read.flags.is_reverse);
//                            image_matrix[base_index][feature_index] += 1;
//                        }
//                    }
                }
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::REF_SKIP:
            case CIGAR_OPERATIONS::PAD:
            case CIGAR_OPERATIONS::DEL:
                // process delete allele here

                if (ref_position -1 >= ref_start && ref_position - 1 <= ref_end) {
                    char ref_base = reference_sequence[ref_position - 1 - ref_start];
                    int base_index = (int)(ref_position - 1 - ref_start + cumulative_observed_insert[ref_position - 1 - ref_start]);
                    int delete_count_index =  get_feature_index('D', read.flags.is_reverse);
                    image_matrix[base_index][delete_count_index] += 1.0;

//                    int feature_index = get_feature_index('*', read.flags.is_reverse);
//                    image_matrix[base_index][feature_index] += 1.0;
                }
                // dont' expand to the full delete length, rather just mount everything to the anchor

                for (int i = 0; i < cigar.length; i++) {
                    if (ref_position + i >= ref_start && ref_position + i <= ref_end) {
                        // update the summary of base
                        int base_index = (int) (ref_position - ref_start + i + cumulative_observed_insert[ref_position - ref_start + i]);
                        char ref_base = reference_sequence[ref_position - ref_start + i];
                        int feature_index = get_feature_index('*', read.flags.is_reverse);

                        image_matrix[base_index][feature_index] += 1.0;
//                        coverage_vector[ref_position - ref_start + i] += 1.0;
                    }
                }

                ref_position += cigar.length;
                break;
            case CIGAR_OPERATIONS::SOFT_CLIP:
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::HARD_CLIP:
                break;
        }
    }
}

RegionalImageSummary RegionalSummaryGenerator::generate_summary(vector <type_read> &reads,
                                                                int chunk_overlap,
                                                                int smaller_chunk_size,
                                                                int feature_size,
                                                                int chunk_id_start,
                                                                bool train_mode) {
    RegionalImageSummary summary{};

    int region_size = (int) (ref_end - ref_start + total_observered_insert_bases + 1);

    // Generate a cover vector of chunk size. Chunk size = 10kb defining the region
    int coverage_vector[ref_end - ref_start + 1];

    // generate the image matrix of chunk_size (10kb) * feature_size (10)
    int** image_matrix = new int*[region_size + 1];

    for (int i = 0; i < region_size + 1; i++) {
        image_matrix[i] = new int[feature_size];
        for (int j = 0; j < feature_size; j++)
            image_matrix[i][j] = 0;
    }

    memset(coverage_vector, 0, sizeof(coverage_vector));

    encode_reference_bases(image_matrix);

    // now iterate over all of the reads and populate the image matrix and coverage matrix
    for (auto &read:reads) {
        // this populates base_summaries and insert_summaries dictionaries
        if(read.mapping_quality > 0) {
            populate_summary_matrix(image_matrix, coverage_vector, read);
        }
    }
    // once the image matrix is generated, scale the counted values.
    for(int i=0;i<region_size;i++){
        // normalize things
//        for(int j=ImageOptionsRegion::BASE_INDEX_START; j < ImageOptionsRegion::BASE_INDEX_START + ImageOptionsRegion::BASE_INDEX_SIZE ; j++){
//            image_matrix[i][j] = (int) (((double)image_matrix[i][j] / max(1.0, (double) coverage_vector[positions[i]-ref_start])) * ImageOptionsRegion::MAX_COLOR_VALUE);
//        }

        for(int j=ImageOptionsRegion::BASE_INDEX_START; j < ImageOptionsRegion::BASE_INDEX_START + ImageOptionsRegion::BASE_INDEX_SIZE ; j++){
            image_matrix[i][j] = (int) min(image_matrix[i][j], ImageOptionsRegion::MAX_COLOR_VALUE);
        }
//        int fwd_feature_index = get_feature_index(ref_at_labels[i], false);
//        int rev_feature_index = get_feature_index(ref_at_labels[i], true);
//
//        for(int j=ImageOptionsRegion::BASE_INDEX_START; j < ImageOptionsRegion::BASE_INDEX_START + ImageOptionsRegion::BASE_INDEX_SIZE ; j++) {
//            if(image_matrix[i][j] == 0) continue;
//
//            if (j == fwd_feature_index || j == rev_feature_index){
//                image_matrix[i][j] = min(ImageOptionsRegion::MAX_REF_COLOR_VALUE, image_matrix[i][j]);
//            } else {
//                image_matrix[i][j] = min(ImageOptionsRegion::MAX_COLOR_VALUE, ImageOptionsRegion::MISMATCH_COLOR_START + image_matrix[i][j]);
//            }
//        }
    }

    labels.resize(region_size + 1, 0);
    labels_variant_type.resize(region_size + 1, 0);
    // check if train mode, if yes, then generate labels
    if(train_mode) {
        for (int i = 0; i < labels_hp1.size(); i++) {
            labels[i] = get_label_index(labels_hp1[i], labels_hp2[i]);
            labels_variant_type[i] = get_variant_type_label_index(variant_type_labels_hp1[i], variant_type_labels_hp2[i]);
        }
    }

    // Now chunk the image into 1kb smaller matrices. Making them 10 * 1000 * feature_size
    int64_t chunk_start = 0;
    int64_t chunk_end = smaller_chunk_size;
    int total_chunks = 0;

    // first count how many chunks will there be
    while(true) {
        total_chunks += 1;
        if (chunk_end == (int) region_size) {
            break;
        }
        chunk_start = chunk_end - chunk_overlap;
        chunk_end = min((int64_t) region_size, chunk_start + smaller_chunk_size);

    }
    // once we know how many chunks there will be, we can generate a pre-defined vector and do the chunking
    summary.chunked_image_matrix.resize(total_chunks, vector<vector<uint8_t> >(smaller_chunk_size,vector<uint8_t>(feature_size)));
    summary.chunked_labels.resize(total_chunks, vector<uint8_t>(smaller_chunk_size));
    summary.chunked_type_labels.resize(total_chunks, vector<uint8_t>(smaller_chunk_size));
    summary.chunked_positions.resize(total_chunks, vector<int64_t>(smaller_chunk_size));
    summary.chunked_index.resize(total_chunks, vector<int32_t>(smaller_chunk_size));
    summary.chunked_ids.resize(total_chunks);

    chunk_start = 0;
    chunk_end = smaller_chunk_size;
    int current_chunk = 0;

    while(true) {
        summary.chunked_ids[current_chunk] = chunk_id_start + current_chunk;
        for(int64_t i = chunk_start; i < chunk_end; i++) {

            if(i<region_size) {
                summary.chunked_positions[current_chunk][i - chunk_start] = positions[i];
                summary.chunked_index[current_chunk][i - chunk_start] = index[i];
                summary.chunked_labels[current_chunk][i - chunk_start] = labels[i];
                summary.chunked_type_labels[current_chunk][i - chunk_start] = labels_variant_type[i];
            } else {
                summary.chunked_positions[current_chunk][i - chunk_start] = -1;
                summary.chunked_index[current_chunk][i - chunk_start] = -1;
                summary.chunked_labels[current_chunk][i - chunk_start] = 0;
                summary.chunked_type_labels[current_chunk][i - chunk_start] = 0;
            }
            for(int j=0; j < feature_size; j++) {
                if(i<=region_size) {
                    summary.chunked_image_matrix[current_chunk][i - chunk_start][j] = image_matrix[i][j];
                } else {
                    summary.chunked_image_matrix[current_chunk][i - chunk_start][j] = 0;
                }
            }
        }
        if (chunk_end >= region_size) {
            break;
        }

        chunk_start = chunk_end - chunk_overlap;
        chunk_end = chunk_start + smaller_chunk_size;

        // in train mode, we don't go outside the window
        if(train_mode) {
            if(chunk_end > region_size) {
                int padding_required = (int)(chunk_end - region_size);
                chunk_start -= padding_required;
                chunk_end = chunk_start + smaller_chunk_size;
                if(chunk_start < 0) break;
            }
        }
        current_chunk += 1;
    }

//    debug_print_matrix(image_matrix, train_mode);
    for (int i = 0; i < region_size + 1; i++)
        free(image_matrix[i]);
    free(image_matrix);
//    exit(0);
//
    return summary;
}


void RegionalSummaryGenerator::debug_print_matrix(int** image_matrix, bool train_mode) {
    cout << "------------- IMAGE MATRIX" << endl;

    cout << setprecision(3);
    for (long long i = ref_start; i <= ref_end; i++) {
        if(i==ref_start) cout<<"REF:\t";
        cout << "  " << reference_sequence[i - ref_start] << "\t";
        if (max_observed_insert[i - ref_start] > 0) {
            for (uint64_t ii = 0; ii < max_observed_insert[i - ref_start]; ii++)cout << "  *" << "\t";
        }
    }
    cout << endl;

    if(train_mode) {
        for (int i = 0; i <= ref_end - ref_start + total_observered_insert_bases; i++) {
            if (i == 0) cout << "TRH1:\t";
            cout << "  " << labels_hp1[i] << "\t";
        }
        cout << endl;

        for (int i = 0; i <= ref_end - ref_start + total_observered_insert_bases; i++) {
            if (i == 0) cout << "TRH2:\t";
            cout << "  " << labels_hp2[i] << "\t";
        }
        cout << endl;
    }

    for (int i = 0; i < labels.size(); i++) {
        if(i==0) cout<<"LBL:\t";
        printf("%3d\t", labels[i]);
    }
    cout << endl;

    for (int i = 0; i < labels_variant_type.size(); i++) {
        if(i==0) cout<<"TYP:\t";
        printf("%3d\t", labels_variant_type[i]);
    }
    cout << endl;

    cout<<"POS:\t";
    for(int i=0; i < positions.size(); i++ ) {
        printf("%3lld\t", positions[i] % 100);
    }
    cout << endl;
    int image_size = ImageOptionsRegion::REFERENCE_INDEX_SIZE + ImageOptionsRegion::BASE_INDEX_SIZE;
    for (int i = 0; i < image_size; i++) {
        cout<< ImageOptionsRegion::column_values[i] <<"\t";
        int region_size = (int) (ref_end - ref_start + total_observered_insert_bases + 1);

        for (int j = 0; j < region_size; j++) {
            printf("%3d\t", image_matrix[j][i]);
        }
        cout << endl;
    }

}