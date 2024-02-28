//
// Created by Kishwar Shafin on 3/21/21.
//

#include "region_summary.h"

#include <utility>
#include "colors.h"
#include <rapidfuzz/fuzz.hpp>

bool fahmid_check = false;

RegionalSummaryGenerator::RegionalSummaryGenerator(string contig, long long region_start, long long region_end, string reference_sequence) {
    this->contig = contig;
    this->ref_start = region_start;
    this->ref_end = region_end;
    this->reference_sequence = std::move(reference_sequence);
    this->total_observered_insert_bases = 0;
    this->max_observed_insert.resize(region_end-region_start+1, 0);
    this->cumulative_observed_insert.resize(region_end-region_start+1, 0);
    fill(this->max_observed_insert.begin(), this->max_observed_insert.end(), 0);
    fill(this->cumulative_observed_insert.begin(), this->cumulative_observed_insert.end(), 0);
    this->print_colored_debug = false;

    if (print_colored_debug){
        cerr << RED << "[REGIONAL SUMMARY] chromosome_name: " << contig << ":" << region_start << "-" << region_end << RESET <<endl;
        cerr << RED << "[REGIONAL SUMMARY] Reference sequence length: " << this->reference_sequence.size() << RESET <<endl;
        cerr << RED << "[REGIONAL SUMMARY] Region range: " <<this->ref_start << "-" << this->ref_end << RESET <<endl;
    }

    this->skip_snps = true;
    this->generate_negatives = false;
    this->candidate_length_thresh = 50;

    this->del_merge_window = 500;
    this->ins_merge_window = 150;
    this->ins_merge_score_threshold = 90;
    this->del_merge_score_threshold = 70;
}

void RegionalSummaryGenerator::generate_max_insert_observed(const type_read& read) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    long long reference_index;
    // if (print_colored_debug){
    //     cerr << MAGENTA << "Read pos: " << read.pos <<" Read ID: " << read.read_id << RESET << endl;
    // }
    if(ImageOptionsRegion::GENERATE_INDELS == true) {
        for (auto &cigar: read.cigar_tuples) {
            cerr << MAGENTA << cigar.operation << " " << cigar.length << RESET << endl;
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
    vector<string> base_labels {"RR", "RA", "RC", "RT", "RG", "R*", "R#", "AA", "AC", "AT", "AG", "A*", "A#", "CC", "CT", "CG", "C*", "C#", "TT", "TG", "T*", "T#", "GG", "G*", "G#", "**", "*#", "##"};


    for(int i=0; i<base_labels.size(); i++) {
        if(base_h1 == base_labels[i][0] && base_h2 == base_labels[i][1]) return i;
        if(base_h2 == base_labels[i][0] && base_h1 == base_labels[i][1]) return i;
    }
    cout<<"TYPE NOT FOUND FOR: "<<base_h1<<" "<<base_h2<<endl;
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

int RegionalSummaryGenerator::get_reference_feature_value(char base) {
    base = toupper(base);
    if (base == 'A') return 1;
    if (base == 'C') return 2;
    if (base == 'G') return 3;
    if (base == 'T') return 4;
    return 5;
}

void RegionalSummaryGenerator::encode_reference_bases(vector< vector<int> >& image_matrix) {
    for (long long ref_position = ref_start; ref_position <= ref_end; ref_position++) {
        // encode the C base
        int base_index = (int) (ref_position - ref_start + cumulative_observed_insert[ref_position - ref_start]);
        // int feature_index = get_reference_feature_index(reference_sequence[ref_position - ref_start]);
        int feature_index = 0;
        int value = get_reference_feature_value(reference_sequence[ref_position - ref_start]);
        image_matrix[base_index][feature_index] = value;

        for(int i = 1; i <= max_observed_insert[ref_position - ref_start]; i++) {
            base_index = (int) (ref_position - ref_start + cumulative_observed_insert[ref_position - ref_start]) + i;
            // feature_index = get_reference_feature_index('*');
             cout<<"\t\t\t\t\t\t\t"<<base_index<<endl;
            feature_index = 0;
            value = get_reference_feature_value(reference_sequence[ref_position - ref_start]);
            image_matrix[base_index][feature_index] = value;
        }
    }
}

bool check_ref_base(char base) {
    if(base=='A' || base=='a' ||
       base=='C' || base=='c' ||
       base=='T' || base=='t' ||
       base =='G' || base=='g') return true;
    return false;
}

int checkIfPreviousSCCandidate(int region_index, const vector< map<string, int> >& AlleleFrequencyMap, const string& candidate_string, bool isDel) {
    // Ensure region_index is within the bounds of AlleleFrequencyMap
    if(region_index < 0 || region_index >= AlleleFrequencyMap.size()) {
        cerr << "Error: region index out of bound" << endl;
        // Handle error: region_index is out of bounds
        return -1;
    }

    if (fahmid_check) {
        cerr << "reg index: " << region_index << endl;
    }

    // Compute left and right limits
    int left_limit = max(0, region_index - 150);
    int right_limit = min(static_cast<int>(AlleleFrequencyMap.size()) - 1, region_index + 150);
    
    // Iterate from left_limit to right_limit inclusive
    for (int i = left_limit; i <= right_limit; ++i) {
        // Check if the candidate_string exists in the current region of AlleleFrequencyMap
        auto it = AlleleFrequencyMap[i].find(candidate_string);
        
        if (it != AlleleFrequencyMap[i].end()) {
            return i;
        }
    }

    if(isDel) {
        for (int i = left_limit; i <= right_limit; ++i) {
            
            string modified_candidate_string = "E" + candidate_string;
            auto it = AlleleFrequencyMap[i].find(modified_candidate_string);
    
            if (it != AlleleFrequencyMap[i].end()) {
                return i;
            }
        }
    }
    return -1;
}

void checkLabels(const vector<char>& labels_hp, int base_index, 
                 const vector<vector<type_truth_record>>& truth_alleles,
                 bool& hasAsterisk, bool& hasHash, string& alt, long& len, int& pos, string& ref) {
    
    int left_limit = max(0, base_index - 150);
    int right_limit = min((int)labels_hp.size() - 1, base_index + 150);
    
    if (fahmid_check) {
        cerr << "left: " << left_limit << endl;
        cerr << "right: " << right_limit << endl;
    }

    hasAsterisk = false;
    hasHash = false;
    alt = "";
    len = 0;
    ref = "";
    for (int i = left_limit; i <= right_limit; ++i) {
        if (labels_hp[i] == '*') {
            hasAsterisk = true;
            int maxLen = -1;
            for (const auto& allele : truth_alleles[i]) {
                long alleleLength = allele.pos_end - allele.pos_start;
                if (alleleLength > maxLen) {
                    maxLen = alleleLength;
                    ref = allele.ref;
                    alt = allele.alt;
                    len = alleleLength;              
                }
            }
            pos = i;
            break;
        } else if (labels_hp[i] == '#') {
            hasHash = true;
            int maxLen = -1;
            // cerr << truth_alleles[i] << endl;
            for (const auto& allele : truth_alleles[i]) {
                long alleleLength = allele.pos_end - allele.pos_start;
                if (alleleLength > maxLen) {
                    maxLen = alleleLength;
                    ref = allele.ref;
                    alt = allele.alt;
                    len = alleleLength;              
                }
            }
            pos = i;
            break;
        }
    }

    if (fahmid_check) {
        cerr << BLUE << "in Check Labels:" << endl;
        cerr << "Length: " << len << endl;
        cerr << "REF: " << ref << endl << RESET;

    }
}


int RegionalSummaryGenerator::get_feature_index(char ref_base, char base, bool is_reverse) {
    base = toupper(base);
    ref_base = toupper(ref_base);
    bool valid_ref = check_ref_base(ref_base);
    if(valid_ref) {
        // this is a mismatch situation
        if (!is_reverse) {
            int start_index  = 7;
            if (base == 'A') return start_index + 1;
            if (base == 'C') return start_index + 2;
            if (base == 'G') return start_index + 3;
            if (base == 'T') return start_index + 4;
            if (base == 'I') return start_index + 5;
            if (base == 'D') return start_index + 6;
            if (base == 'S') return 26;
            return start_index + 7;
        } else {
            // tagged and forward
            int start_index  = 18;
            if (base == 'A') return start_index + 1;
            if (base == 'C') return start_index + 2;
            if (base == 'G') return start_index + 3;
            if (base == 'T') return start_index + 4;
            if (base == 'I') return start_index + 5;
            if (base == 'D') return start_index + 6;
            if (base == 'S') return 27;
            return start_index + 7;
        }
    } else {
        return -1;
    }
}



void RegionalSummaryGenerator::generate_labels(const vector<type_truth_record>& hap1_records, const vector<type_truth_record>& hap2_records) {
    int region_size = (int) (ref_end - ref_start + total_observered_insert_bases + 1);
    labels_hp1.resize(region_size + 1, '*');
    labels_hp2.resize(region_size + 1, '*');
    variant_type_labels_hp1.resize(region_size + 1, VariantTypes::HOM_REF);
    variant_type_labels_hp2.resize(region_size + 1, VariantTypes::HOM_REF);
    // Fixed labels for structural variants

    hp1_truth_alleles.resize(region_size + 1);
    hp2_truth_alleles.resize(region_size + 1);

    for(long long pos = ref_start; pos <= ref_end; pos++) {
        int base_index = (int)(pos - ref_start + cumulative_observed_insert[pos - ref_start]);
        labels_hp1[base_index] = 'R';
        labels_hp2[base_index] = 'R';
    }

    for(const auto& truth_record : hap1_records) {
        if (print_colored_debug){
            cerr << BLUE <<"TRUTH: "<<truth_record.contig <<" "<<truth_record.pos_start << "-" << truth_record.pos_end <<" "<<truth_record.ref <<"--" << truth_record.ref.size() << " " <<truth_record.alt<< RESET <<endl;
        }
        if(truth_record.ref.length() > truth_record.alt.length()) {
            //it's a delete
           //cerr << BLUE <<"TRUTH DELETE: "<<truth_record.contig<<" "<<truth_record.pos_start<<" "<<truth_record.ref<<" "<<truth_record.alt<< RESET <<endl;
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp1[base_index] = VariantTypes::DELETE;
                labels_hp1[base_index] = '#';
                hp1_truth_alleles[base_index].push_back(truth_record);

                long long cum_ob_ins_vect_idx = min(truth_record.pos_end - ref_start - 1, (long long)cumulative_observed_insert.size() - 1);
                int base_end_index = (int) (truth_record.pos_end - ref_start - 1 + cumulative_observed_insert[cum_ob_ins_vect_idx]);
                base_end_index = min(base_end_index, region_size);
                variant_type_labels_hp1[base_end_index] = VariantTypes::DELETE;
                labels_hp1[base_end_index]='#';
                hp1_truth_alleles[base_end_index].push_back(truth_record);

                if(print_colored_debug){
                    cerr << "labels_hp1: " << base_index << " " << labels_hp1[base_index] << endl;
                    cerr << base_end_index << " " << labels_hp1[base_end_index] << endl;
                }
            }
        } else if(truth_record.ref.length() < truth_record.alt.length()) {
            //it's an insert
//            cout<<"TRUTH INSERT: "<<truth_record.contig<<" "<<truth_record.pos_start<<" "<<truth_record.ref<<" "<<truth_record.alt<<endl;
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp1[base_index] = VariantTypes::INSERT;
                labels_hp1[base_index] = '*';
                hp1_truth_alleles[base_index].push_back(truth_record);
            }
        } else if(truth_record.ref.length() == truth_record.alt.length()) {
            //it's a SNP
        //    cerr<<"TRUTH SNP: "<<truth_record.contig<<" "<<truth_record.pos_start<<" "<<truth_record.ref<<" "<<truth_record.alt<<endl;
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp1[base_index] = VariantTypes::SNP;
                hp1_truth_alleles[base_index].push_back(truth_record);
            }

            for(long long pos = truth_record.pos_start; pos < truth_record.pos_end; pos++) {
                if (pos >= ref_start && pos <= ref_end) {
                    int base_index = (int) (pos - ref_start + cumulative_observed_insert[pos - ref_start]);
                    char ref_base = reference_sequence[pos - ref_start];
                    char alt_base = truth_record.alt[pos - truth_record.pos_start];
                    if(ref_base == alt_base) {
                        labels_hp1[base_index] = 'R';
                    } else {
                        labels_hp1[base_index] = truth_record.alt[pos - truth_record.pos_start];
                        // cerr << "labels_hp1: " << base_index << " " << labels_hp1[base_index] << endl;
                    }
                }
            }
        }
    }

    for(const auto& truth_record : hap2_records) {
        if(truth_record.ref.length() > truth_record.alt.length()) {
            if (print_colored_debug){
                cerr << BLUE <<"TRUTH: "<<truth_record.contig <<" "<<truth_record.pos_start << "-" << truth_record.pos_end <<" "<<truth_record.ref <<"--" << truth_record.ref.size() << " " <<truth_record.alt<< RESET <<endl;
            }
            //it's a delete
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp2[base_index] = VariantTypes::DELETE;
                labels_hp2[base_index] = '#';
                hp2_truth_alleles[base_index].push_back(truth_record);

                long long cum_ob_ins_vect_idx = min(truth_record.pos_end - ref_start - 1, (long long)cumulative_observed_insert.size() - 1);
                int base_end_index = (int) (truth_record.pos_end - ref_start - 1 + cumulative_observed_insert[cum_ob_ins_vect_idx]);
                base_end_index = min(base_end_index, region_size);
                variant_type_labels_hp2[base_end_index] = VariantTypes::DELETE;
                labels_hp2[base_end_index]='#';
                hp2_truth_alleles[base_end_index].push_back(truth_record);
                if (print_colored_debug){
                    cerr << "labels_hp2: " << base_index << " " << labels_hp2[base_index] << endl;
                    cerr << base_end_index << " " << labels_hp2[base_end_index] << endl;
                }
            }
        } else if(truth_record.ref.length() < truth_record.alt.length()) {
            //it's an insert
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp2[base_index] = VariantTypes::INSERT;
                labels_hp2[base_index] = '*';
                hp2_truth_alleles[base_index].push_back(truth_record);
            }
        } else if(truth_record.ref.length() == truth_record.alt.length()) {
            //it's a SNP
            // cerr<<"TRUTH SNP: "<<truth_record.contig<<" "<<truth_record.pos_start<<" "<<truth_record.ref<<" "<<truth_record.alt<<endl;
            if (truth_record.pos_start >= ref_start && truth_record.pos_start <= ref_end) {
                int base_index = (int) (truth_record.pos_start - ref_start + cumulative_observed_insert[truth_record.pos_start - ref_start]);
                variant_type_labels_hp2[base_index] = VariantTypes::SNP;
                hp2_truth_alleles[base_index].push_back(truth_record);
            }

            for(long long pos = truth_record.pos_start; pos < truth_record.pos_end; pos++) {
                if (pos >= ref_start && pos <= ref_end) {
                    int base_index = (int) (pos - ref_start + cumulative_observed_insert[pos - ref_start]);

                    char ref_base = reference_sequence[pos - ref_start];
                    char alt_base = truth_record.alt[pos - truth_record.pos_start];
                    if(ref_base == alt_base) {
                        labels_hp2[base_index] = 'R';
                    } else {
                        labels_hp2[base_index] = truth_record.alt[pos - truth_record.pos_start];
                    }
                }
            }
        }
    }
}


void RegionalSummaryGenerator::populate_summary_matrix(vector< vector<int> >& image_matrix,
                                                       int *coverage_vector,
                                                       int *snp_count,
                                                       int *insert_count,
                                                       int *delete_count,
                                                       vector< map<string, int> > &AlleleFrequencyMap,
                                                       vector< map<string, int> > &AlleleFrequencyMapFwdStrand,
                                                       vector< map<string, int> > &AlleleFrequencyMapRevStrand,
                                                       vector< set<string> > &AlleleMap,
                                                       type_read read,
                                                       double min_snp_baseq,
                                                       double min_indel_baseq,
                                                       bool train_mode) {

    // cerr << "Populate Summary Matrix Called" << endl;
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    double base_quality = read.base_qualities[read_index];
    // bool debug_populate_summary = (read.query_name == "SRR9001772.257834");
    bool debug_populate_summary = (read.query_name == "SRR9001772.906641");
    char cigarOpMap[11] = {
        [CIGAR_OPERATIONS::MATCH] = 'M',
        [CIGAR_OPERATIONS::IN] = 'I',
        [CIGAR_OPERATIONS::DEL] = 'D',
        [CIGAR_OPERATIONS::REF_SKIP] = 'N',
        [CIGAR_OPERATIONS::SOFT_CLIP] = 'S',
        [CIGAR_OPERATIONS::HARD_CLIP] = 'H',
        [CIGAR_OPERATIONS::PAD] = 'P',
        [CIGAR_OPERATIONS::EQUAL] = '=',
        [CIGAR_OPERATIONS::DIFF] = 'X',
        [CIGAR_OPERATIONS::BACK] = 'B',
        [10] = '!'
    };
    for (int cigar_i=0; cigar_i<read.cigar_tuples.size(); cigar_i++) {
        CigarOp cigar = read.cigar_tuples[cigar_i];

        if (ref_position > ref_end) {
            if (print_colored_debug){
                cerr << BOLDRED << "breaking" << RESET << endl;
                cerr << "Ref position: " << ref_position << endl;
            }
            break;
        }
        char cop = cigarOpMap[cigar.operation==CIGAR_OPERATIONS::UNSPECIFIED?10:cigar.operation];

        if (print_colored_debug){
            // cerr << BOLDWHITE << "cigar: " << cop << " " << cigar.length << RESET << endl;
        }
        switch (cigar.operation) {
            case CIGAR_OPERATIONS::EQUAL:
            case CIGAR_OPERATIONS::DIFF:
            case CIGAR_OPERATIONS::MATCH:
                // cerr << "On Match" << endl;
                cigar_index = 0;
                if (ref_position < ref_start) {
                    cigar_index = min(ref_start - ref_position, (long long) cigar.length);
                    read_index += cigar_index;
                    ref_position += cigar_index;
                }
                for (int i = cigar_index; i < cigar.length; i++) {
                    base_quality = read.base_qualities[read_index];

                    if (ref_position >= ref_start && ref_position <= ref_end) {
                        char base = read.sequence[read_index];
                        char ref_base = reference_sequence[ref_position - ref_start];
                        string alt(1, read.sequence[read_index]);

                        int base_index = (int)(ref_position - ref_start + cumulative_observed_insert[ref_position - ref_start]);
                        int feature_index = get_feature_index(ref_base, base, read.flags.is_reverse);

                        // update the summary of base
                        if(base_quality >= min_snp_baseq) {
                            coverage_vector[ref_position - ref_start] += 1;
                            // look front and see if it's anchoring an INSERT or DELETE
                            if(i == cigar.length - 1 && cigar_i != read.cigar_tuples.size() - 1) {
                                CigarOp next_cigar = read.cigar_tuples[cigar_i + 1];
                                if(next_cigar.operation != CIGAR_OPERATIONS::IN && next_cigar.operation != CIGAR_OPERATIONS::DEL) {
                                    if(!read.flags.is_reverse) image_matrix[base_index][4] -= 1;
                                    else image_matrix[base_index][15] -= 1;
                                }
                            }
                            else {
                                if (!read.flags.is_reverse) image_matrix[base_index][4] -= 1;
                                else image_matrix[base_index][15] -= 1;
                            }
                        }

                        if(ref_base != base && base_quality >= min_snp_baseq) {
                        // For picking negative examples
                        // if (ref_position == (ref_start + ref_end) / 2){
                            snp_count[ref_position - ref_start] += 1;
                            if(feature_index >= 0) image_matrix[base_index][feature_index] -= 1;
                            // save the candidate
                            string candidate_string = char(AlleleType::SNP_ALLELE + '0') + alt;

                            int region_index = (int) (ref_position - ref_start);

                            if (AlleleFrequencyMap[region_index].find(candidate_string) != AlleleFrequencyMap[region_index].end()) {
                                AlleleFrequencyMap[region_index][candidate_string] += 1;
                                if(read.flags.is_reverse) {
                                    AlleleFrequencyMapRevStrand[region_index][candidate_string] += 1;
                                } else {
                                    AlleleFrequencyMapFwdStrand[region_index][candidate_string] += 1;
                                }
                            } else {
                                AlleleFrequencyMap[region_index][candidate_string] = 1;
                                if(read.flags.is_reverse) {
                                    AlleleFrequencyMapFwdStrand[region_index][candidate_string] = 0;
                                    AlleleFrequencyMapRevStrand[region_index][candidate_string] = 1;
                                } else {
                                    AlleleFrequencyMapFwdStrand[region_index][candidate_string] = 1;
                                    AlleleFrequencyMapRevStrand[region_index][candidate_string] = 0;
                                }
                            }

                            if (AlleleMap[region_index].find(candidate_string) == AlleleMap[region_index].end())
                                AlleleMap[region_index].insert(candidate_string);
                        } else if (base_quality >= min_snp_baseq){
                            if(feature_index >= 0)  image_matrix[base_index][feature_index] -= 1;

                        }
                    }
                    read_index += 1;
                    ref_position += 1;
                }
                break;
            case CIGAR_OPERATIONS::IN:
                // cerr << "On Insertion" << endl;
                if (ref_position - 1 >= ref_start && ref_position - 1 <= ref_end && read_index - 1 >= 0) {
                    // process insert allele here
                    string alt;
                    char ref_base = reference_sequence[ref_position - 1 - ref_start];
                    int base_index = (int)((ref_position - 1) - ref_start + cumulative_observed_insert[(ref_position - 1) - ref_start]);
                    int insert_count_index =  get_feature_index(ref_base, 'I', read.flags.is_reverse);

                    if (read_index - 1 >= 0) alt = read.sequence.substr(read_index - 1, cigar.length + 1);
                    else alt = ref_base + read.sequence.substr(read_index, cigar.length);

                    int len = cigar.length + 1;
                    base_quality = 0;
                    int start_index = 0;
                    if (read_index - 1 >= 0) start_index  = read_index - 1;
                    else start_index = read_index;

                    for(int i = start_index; i < start_index + len; i++) {
                        base_quality += read.base_qualities[i];
                    }

                    // include reads that were excluded due to anchor bases' quality
                    if(base_quality >= min_indel_baseq * len && read.base_qualities[start_index] < min_snp_baseq)
                        coverage_vector[ref_position - 1 - ref_start] += 1;


                    // save the candidate
                    string candidate_string = char(AlleleType::INSERT_ALLELE + '0') + alt;

                    // only process candidates that are smaller than 50bp as they 50bp+ means SV
                    if(candidate_string.length() >= candidate_length_thresh && base_quality >= min_indel_baseq * len) {
                        // if(insert_count_index >= 0) image_matrix[base_index][insert_count_index] -= 1;
                        // insert_count[ref_position - 1 - ref_start] += 1;
                        int region_index = (int) (ref_position - 1 - ref_start);

                        // Outer for loop scans the region for 16 bases on either side of the candidate
                        // This is to merge the scattered candidates that are close to each other into a single candidate
                        for(int i = max(0, region_index - ins_merge_window); i < min(int(AlleleMap.size()), region_index + ins_merge_window); i++){
                            bool found = false;
                            // iterate over AlleleMap[region_index] and check if there is a match
                            for(auto it = AlleleMap[i].begin(); it != AlleleMap[i].end(); it++){
                                string existing_candidate = *it;
                                double score = rapidfuzz::fuzz::ratio(existing_candidate, candidate_string);
                                if (score >= ins_merge_score_threshold && existing_candidate[0] == candidate_string[0]){
                                    candidate_string = existing_candidate;
                                    region_index = i;
                                    base_index = (int) (i+cumulative_observed_insert[i]);
                                    found = true;
                                    if (print_colored_debug)
                                    cerr << "Found potential match INS: " << score << " " << region_index + ref_start + 1 << endl;
                                    break;
                                }
                            }
                            if (found) break;
                        }
                        insert_count[region_index] += 1;
                        if(insert_count_index >= 0) image_matrix[base_index][insert_count_index] -= 1;

                        if (AlleleFrequencyMap[region_index].find(candidate_string) != AlleleFrequencyMap[region_index].end()) {
                            AlleleFrequencyMap[region_index][candidate_string] += 1;
                            if(read.flags.is_reverse) {
                                AlleleFrequencyMapRevStrand[region_index][candidate_string] += 1;
                            } else {
                                AlleleFrequencyMapFwdStrand[region_index][candidate_string] += 1;
                            }
                        } else {
                            AlleleFrequencyMap[region_index][candidate_string] = 1;
                            if(read.flags.is_reverse) {
                                AlleleFrequencyMapFwdStrand[region_index][candidate_string] = 0;
                                AlleleFrequencyMapRevStrand[region_index][candidate_string] = 1;
                            } else {
                                AlleleFrequencyMapFwdStrand[region_index][candidate_string] = 1;
                                AlleleFrequencyMapRevStrand[region_index][candidate_string] = 0;
                            }
                        }

                        if (AlleleMap[region_index].find(candidate_string) == AlleleMap[region_index].end())
                            AlleleMap[region_index].insert(candidate_string);
                    }
                }
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::DEL:
                // cerr << "On Delete" << endl;
                // process delete allele here
                if (ref_position -1 >= ref_start && ref_position - 1 <= ref_end) {
                    int region_index = (int) (ref_position - 1 - ref_start);
                    int region_end_index = (int) (ref_position - ref_start + cigar.length - 1);
                    char ref_base = reference_sequence[region_index];
                    int base_index = (int)(region_index + cumulative_observed_insert[region_index]);
                    char ref_end_base = reference_sequence[region_end_index];
                    int base_end_index = (int)(region_end_index + cumulative_observed_insert[region_end_index]);
                
                    int delete_count_index =  get_feature_index(ref_base, 'D', read.flags.is_reverse);
                    int delete_count_end_index = get_feature_index(ref_end_base, 'D', read.flags.is_reverse);
                    
                    // Helpful debug information
                    if (print_colored_debug) {
                        // cerr << HCYN << "ref_index: " << ref_position - 1 - ref_start << " region_end_index: " << region_end_index << endl;
                        // cerr << "base_index: " << base_index << " base_end_index: " << base_end_index << endl;
                        // cerr << "ref_base: " << ref_base << " ref_end_base: " << ref_end_base << endl;
                        // cerr << "del_count_idx: " << delete_count_index << " del_count_end_idx: " << delete_count_end_index << RESET <<endl;
                    }
                    // process delete allele here
                    string ref = reference_sequence.substr(ref_position - ref_start - 1, cigar.length + 1);
                    if (debug_populate_summary){
                        cerr << "DEL: " << ref << endl;
                        cerr << "read_index: " << read_index << endl;
                    }
                    string alt;

                    if (read_index - 1 >= 0) alt = read.sequence.substr(read_index - 1, 1);
                    else alt = reference_sequence.substr(ref_position - ref_start - 1, 1);

                    base_quality = read.base_qualities[read_index];
                    string candidate_string = char(AlleleType::DELETE_ALLELE + '0') + ref;

                    // only process candidates that are smaller than 50bp as they 50bp+ means SV
                    // no base-quality check for deletes
                    // modified to pick up deletion candidates of larger lengths
                    if(candidate_string.length() >= candidate_length_thresh) {
                        for(int i = max(0, region_index - del_merge_window); i < min(int(AlleleMap.size()), region_index + del_merge_window); i++){
                            bool found = false;
                            // iterate over AlleleMap[region_index] and check if there is a match
                            for(auto it = AlleleMap[i].begin(); it != AlleleMap[i].end(); it++){
                                string existing_candidate = *it;
                                double score = rapidfuzz::fuzz::ratio(existing_candidate, candidate_string);
                                if (score >= del_merge_score_threshold && existing_candidate[0] == candidate_string[0]){
                                    candidate_string = existing_candidate;
                                    if(print_colored_debug) cerr << "Region index: " << ref_position << endl;
                                    region_index = i;
                                    region_end_index = i + candidate_string.length() - 2; // 1 for the first character and 1 for 0-based index
                                    base_index = (int) (i+cumulative_observed_insert[i]);
                                    base_end_index = (int) (region_end_index + cumulative_observed_insert[region_end_index]);
                                    found = true;
                                    if(print_colored_debug)
                                    cerr << "Found potential match DEL: " << score << " " << region_index + ref_start + 1 << " Ending at " << region_end_index + ref_start  << endl;
                                    break;
                                }
                            }
                            if (found) break;
                        }
                        delete_count[region_index] += 1;
                        delete_count[region_end_index] += 1;
                        if(delete_count_index >= 0) image_matrix[base_index][delete_count_index] -= 1.0;
                        if(delete_count_end_index >= 0) image_matrix[base_end_index][delete_count_end_index] -= 1.0;
                        // also increment delete_count at the end of the deletion
                        if (print_colored_debug) {
                            cerr << BOLDWHITE <<  "Deletion end base: " << region_end_index << RESET <<endl;
                        }

                        if (AlleleFrequencyMap[region_index].find(candidate_string) != AlleleFrequencyMap[region_index].end()) {
                            // TODO: Merge candidate strings here
                            AlleleFrequencyMap[region_index][candidate_string] += 1;
                            AlleleFrequencyMap[region_end_index]['E' + candidate_string] += 1;
                            if(read.flags.is_reverse) {
                                AlleleFrequencyMapRevStrand[region_index][candidate_string] += 1;
                                AlleleFrequencyMapRevStrand[region_end_index]['E' + candidate_string] += 1;
                            } else {
                                AlleleFrequencyMapFwdStrand[region_index][candidate_string] += 1;
                                AlleleFrequencyMapFwdStrand[region_end_index]['E' + candidate_string] += 1;
                            }
                        } else {
                            AlleleFrequencyMap[region_index][candidate_string] = 1;
                            AlleleFrequencyMap[region_end_index]['E' + candidate_string] = 1;
                            if(read.flags.is_reverse) {
                                AlleleFrequencyMapFwdStrand[region_index][candidate_string] = 0;
                                AlleleFrequencyMapRevStrand[region_index][candidate_string] = 1;
                                AlleleFrequencyMapFwdStrand[region_end_index]['E' + candidate_string] = 0;
                                AlleleFrequencyMapRevStrand[region_end_index]['E' + candidate_string] = 1;
                            } else {
                                AlleleFrequencyMapFwdStrand[region_index][candidate_string] = 1;
                                AlleleFrequencyMapRevStrand[region_index][candidate_string] = 0;
                                AlleleFrequencyMapFwdStrand[region_end_index]['E' + candidate_string] = 1;
                                AlleleFrequencyMapRevStrand[region_end_index]['E' + candidate_string] = 0;
                            }
                        }

                        if (AlleleMap[region_index].find(candidate_string) == AlleleMap[region_index].end()) {
                            AlleleMap[region_index].insert(candidate_string);
                            AlleleMap[region_end_index].insert('E' + candidate_string);
                        }
                    }
                    // cout<<"DEL: "<<ref_position<<" "<<ref<<" "<<alt<<" "<<AlleleFrequencyMap[candidate_alt]<<endl;
                    for (int i = 0; i < cigar.length; i++) {
                        if ((region_index + 1 + i) >= ref_start && (region_index + 1 + i)<= ref_end) {
                            // update the summary of base
                            int base_index = (int) (region_index + 1 + i + cumulative_observed_insert[region_index + 1 + i]);
                            char ref_base = reference_sequence[region_index + 1 + i];
                            int feature_index = get_feature_index(ref_base, '*', read.flags.is_reverse);

                            if(feature_index >= 0)  image_matrix[base_index][feature_index] -= 1;
    //                        coverage_vector[ref_position - ref_start + i] += 1.0;
                        }
                    }


                }
                // dont' expand to the full delete length, rather just mount everything to the anchor

                ref_position += cigar.length;
                break;
            case CIGAR_OPERATIONS::REF_SKIP:
            case CIGAR_OPERATIONS::PAD:
                ref_position += cigar.length;
            
            case CIGAR_OPERATIONS::SOFT_CLIP:
                if (fahmid_check) {
                    cerr << RED << "Fahmid Check for ref pos detection. REFPOS: " << ref_position << endl <<RESET;
                }

                if (ref_position -1 >= ref_start && ref_position - 1 <= ref_end) {
                    if (fahmid_check) {
                        cerr << GREEN << "Fahmid Check for ref pos detection. REFPOS: " << ref_position << endl <<RESET;
                        cerr << "Fahmid Check ref pos: " << ref_position << endl;
                    }
                    char ref_base = reference_sequence[ref_position - 1 - ref_start];
                    int base_index = (int)(ref_position - 1 - ref_start + cumulative_observed_insert[ref_position - 1 - ref_start]);
                    
                    if (train_mode) {
                        bool hp1_hasAsterisk, hp1_hasHash, hp2_hasAsterisk, hp2_hasHash;
                        string hp1_alt, hp2_alt, hp1_ref, hp2_ref;
                        long len_1, len_2;
                        int pos_1, pos_2;

                        checkLabels(labels_hp1, base_index, hp1_truth_alleles, hp1_hasAsterisk, hp1_hasHash, hp1_alt, len_1, pos_1, hp1_ref);
                        checkLabels(labels_hp2, base_index, hp2_truth_alleles, hp2_hasAsterisk, hp2_hasHash, hp2_alt, len_2, pos_2, hp2_ref);

                        if(hp1_hasAsterisk || hp2_hasAsterisk) {
                            // Insertion Case
                            if (fahmid_check) {
                                cerr << "Found Insertion" <<endl;
                            }

                            string alt;
                            int len, pos;
                            
                            if(hp1_hasAsterisk) {
                                alt = hp1_alt;
                                len = len_1;
                            } else {
                                alt = hp2_alt;
                                len = len_2;
                            }

                            if (fahmid_check) {
                                cerr << len << endl;
                                cerr << alt << endl;
                            }
                            
                            // save the candidate
                            string candidate_string = char(AlleleType::INSERT_ALLELE + '0') + alt;
                            
                            if (fahmid_check) {
                                cerr << "Before reg index func call" << endl;
                            }

                            int region_index = checkIfPreviousSCCandidate(ref_position - 1 - ref_start, AlleleFrequencyMap, candidate_string, false);
                            
                            
                            if(candidate_string.length() >= candidate_length_thresh) {
                                
                                if (region_index == -1) {
                                    region_index = (int) (ref_position - 1 - ref_start);
                                }
                                
                                insert_count[ref_position - 1 - ref_start] += 1;

                                if (AlleleFrequencyMap[region_index].find(candidate_string) != AlleleFrequencyMap[region_index].end()) {
                                    AlleleFrequencyMap[region_index][candidate_string] += 1;
                                    if(read.flags.is_reverse) {
                                        AlleleFrequencyMapRevStrand[region_index][candidate_string] += 1;
                                    } else {
                                        AlleleFrequencyMapFwdStrand[region_index][candidate_string] += 1;
                                    }
                                } else {
                                    AlleleFrequencyMap[region_index][candidate_string] = 1;
                                    if(read.flags.is_reverse) {
                                        AlleleFrequencyMapFwdStrand[region_index][candidate_string] = 0;
                                        AlleleFrequencyMapRevStrand[region_index][candidate_string] = 1;
                                    } else {
                                        AlleleFrequencyMapFwdStrand[region_index][candidate_string] = 1;
                                        AlleleFrequencyMapRevStrand[region_index][candidate_string] = 0;
                                    }
                                }

                                if (AlleleMap[region_index].find(candidate_string) == AlleleMap[region_index].end())
                                    AlleleMap[region_index].insert(candidate_string);
                            }

                        } else if(hp1_hasHash || hp2_hasHash) {
                            // Deletion Case
                            if (fahmid_check) {
                                cerr << "Found Deletion" <<endl;
                            }

                            string ref;
                            int len, pos;
                            if (fahmid_check) {
                                // cerr << "Fahmid Check" << endl;
                            }

                            if(hp1_hasHash) {
                                if (fahmid_check) {
                                    cerr << "at hp1 hash" << endl;
                                }
                                ref = hp1_ref;
                                len = len_1;
                            } else {
                                if (fahmid_check) {
                                    cerr << "at hp2 hash" << endl;
                                }
                                ref = hp2_ref;
                                len = len_2;
                            }
                            if (fahmid_check) {
                                cerr << len << endl;
                                cerr << ref << endl;
                            }

                            string candidate_string = char(AlleleType::DELETE_ALLELE + '0') + ref;

                            int region_index = checkIfPreviousSCCandidate(ref_position - 1 - ref_start, AlleleFrequencyMap, candidate_string, true);
                            if (region_index == -1) {
                                    region_index = (int) (ref_position - 1 - ref_start);
                            }

                            int region_end_index = (int) (region_index + len);
                            
                            int ref_ahead_index = region_index + len;
                            int ref_behind_index = region_index - len;
                            

                            if (fahmid_check) {
                                // cerr << "Fahmid Check" << endl;
                            }

                            if (ref_ahead_index >= AlleleFrequencyMap.size() && ref_behind_index >= 0) {
                                if (fahmid_check) {
                                    cerr << "Fahmid Check: may be left clip" << endl;
                                }
                                region_end_index = region_index;
                                region_index = (int)ref_behind_index;
                                
                            } else if (ref_ahead_index < AlleleFrequencyMap.size() && ref_behind_index < 0) {
                                if (fahmid_check) {
                                    cerr << "Fahmid Check: may be right clip" << endl;
                                }
                            } else {
                                if (fahmid_check) {
                                    cerr << RED << "[" << ref_behind_index << ", " << ref_ahead_index << "]" << AlleleFrequencyMap.size() << endl;
                                    cerr << "Fahmid Check: Some Clip Error!!" << endl << RESET;
                                }
                            }

                            if(candidate_string.length() >= candidate_length_thresh) {
                                if (fahmid_check) {
                                    cerr << "Fahmid Check refpos: " << ref_position << ", refstart: " << ref_start << endl;
                                
                                    // cerr << region_index << endl;
                                    // cerr << region_end_index << endl;

                                    cerr << "Fahmid Check b4 delcount: [" << ref_position - 1 - ref_start << ", " << ref_position - ref_start + len - 1 << "]" << endl;
                                }

                                delete_count[region_index] += 1;

                                delete_count[region_end_index] += 1;

                                if (fahmid_check) {
                                    cerr << "Fahmid Check after delcount" <<endl;
                                }
                                
                                if (AlleleFrequencyMap[region_index].find(candidate_string) != AlleleFrequencyMap[region_index].end()) {
                                    if (fahmid_check) {
                                        cerr << "Fahmid Check in AFM not first" << endl;
                                    }

                                    AlleleFrequencyMap[region_index][candidate_string] += 1;
                                    AlleleFrequencyMap[region_end_index]['E' + candidate_string] += 1;
                                    if(read.flags.is_reverse) {
                                        AlleleFrequencyMapRevStrand[region_index][candidate_string] += 1;
                                        AlleleFrequencyMapRevStrand[region_end_index]['E' + candidate_string] += 1;
                                    } else {
                                        AlleleFrequencyMapFwdStrand[region_index][candidate_string] += 1;
                                        AlleleFrequencyMapFwdStrand[region_end_index]['E' + candidate_string] += 1;
                                    }
                                } else {
                                    if (fahmid_check) {
                                        cerr << "Fahmid Check in AFM first" << endl;
                                    }
                                    AlleleFrequencyMap[region_index][candidate_string] = 1;
                                    if (fahmid_check) {
                                        // cerr << "Fahmid Check " << AlleleFrequencyMap.size() << endl;
                                    }
                                    AlleleFrequencyMap[region_end_index]['E' + candidate_string] = 1;
                                    if(read.flags.is_reverse) {
                                        AlleleFrequencyMapFwdStrand[region_index][candidate_string] = 0;
                                        AlleleFrequencyMapRevStrand[region_index][candidate_string] = 1;
                                        AlleleFrequencyMapFwdStrand[region_end_index]['E' + candidate_string] = 0;
                                        AlleleFrequencyMapRevStrand[region_end_index]['E' + candidate_string] = 1;
                                    } else {
                                        AlleleFrequencyMapFwdStrand[region_index][candidate_string] = 1;
                                        AlleleFrequencyMapRevStrand[region_index][candidate_string] = 0;
                                        AlleleFrequencyMapFwdStrand[region_end_index]['E' + candidate_string] = 1;
                                        AlleleFrequencyMapRevStrand[region_end_index]['E' + candidate_string] = 0;
                                    }
                                }
                                
                                if (fahmid_check) {
                                    // cerr << "Fahmid Check" << endl;
                                }
                                if (AlleleMap[region_index].find(candidate_string) == AlleleMap[region_index].end()) {
                                    AlleleMap[region_index].insert(candidate_string);
                                    AlleleMap[region_end_index].insert('E' + candidate_string);
                                }
                
                                for (int i = 0; i < len; i++) {
                                    if (ref_position + i >= ref_start && ref_position + i <= ref_end) {
                                        // update the summary of base
                                        int base_index = (int) (ref_position - ref_start + i + cumulative_observed_insert[ref_position - ref_start + i]);
                                        char ref_base = reference_sequence[ref_position - ref_start + i];
                                        int feature_index = get_feature_index(ref_base, '*', read.flags.is_reverse);

                                        if(feature_index >= 0)  image_matrix[base_index][feature_index] -= 1;

                                    }
                                }
                            }    
                        
                        }
                    
                    } else {
                        
                        // string ref_base_str{ref_base};
                        // clip cigar length to 100
                        string alt = ref_base + read.sequence.substr(read_index, min(100, cigar.length));
                        string candidate_string = char(AlleleType::INSERT_ALLELE + '0') + alt;
                            
                        if(true) {
                            
                            int region_index = (int) (ref_position - 1 - ref_start);
                            
                            insert_count[ref_position - 1 - ref_start] += 1;

                            if (AlleleFrequencyMap[region_index].find(candidate_string) != AlleleFrequencyMap[region_index].end()) {
                                AlleleFrequencyMap[region_index][candidate_string] += 1;
                                if(read.flags.is_reverse) {
                                    AlleleFrequencyMapRevStrand[region_index][candidate_string] += 1;
                                } else {
                                    AlleleFrequencyMapFwdStrand[region_index][candidate_string] += 1;
                                }
                            } else {
                                AlleleFrequencyMap[region_index][candidate_string] = 1;
                                if(read.flags.is_reverse) {
                                    AlleleFrequencyMapFwdStrand[region_index][candidate_string] = 0;
                                    AlleleFrequencyMapRevStrand[region_index][candidate_string] = 1;
                                } else {
                                    AlleleFrequencyMapFwdStrand[region_index][candidate_string] = 1;
                                    AlleleFrequencyMapRevStrand[region_index][candidate_string] = 0;
                                }
                            }

                            if (AlleleMap[region_index].find(candidate_string) == AlleleMap[region_index].end())
                                AlleleMap[region_index].insert(candidate_string);
                        }

                    }

                    int soft_clip_count_index =  get_feature_index(ref_base, 'S', read.flags.is_reverse);

                    if (fahmid_check) {
                        cerr << "SC INDEX: " << soft_clip_count_index << endl;
                    }

                    if(soft_clip_count_index >= 0) image_matrix[base_index][soft_clip_count_index] -= 1.0;
                    
                    if (fahmid_check){
                        cerr << "SC read_index: " << read_index << endl;
                    }

                }

                read_index += cigar.length;
                break;
            
            case CIGAR_OPERATIONS::HARD_CLIP:
                break;
        }
    }
}

vector<CandidateImageSummary> RegionalSummaryGenerator::generate_summary(vector <type_read> &reads,
                                                                         double min_snp_baseq,
                                                                         double min_indel_baseq,
                                                                         double snp_freq_threshold,
                                                                         double insert_freq_threshold,
                                                                         double delete_freq_threshold,
                                                                         double min_coverage_threshold,
                                                                         double snp_candidate_freq_threshold,
                                                                         double indel_candidate_freq_threshold,
                                                                         double candidate_support_threshold,
                                                                         bool skip_indels,
                                                                         long long candidate_region_start,
                                                                         long long candidate_region_end,
                                                                         int candidate_window_size,
                                                                         int feature_size,
                                                                         bool train_mode) {
    int region_size = (int) (ref_end - ref_start + total_observered_insert_bases + 1);

    if (fahmid_check) {
        cerr << "Fahmid Check Region Size: " << region_size << endl;
    }

    // cerr << "Generate Summary Called" << endl;
    if(print_colored_debug){
        cerr << BLUE << "ref_start, ref_end: [" << ref_start << "-" << ref_end << "] total_observed_insert_bases: " << total_observered_insert_bases << " region_size: " << region_size << RESET << endl;
    }
    // Generate a cover vector of chunk size. Chunk size = 10kb defining the region
    int coverage_vector[ref_end - ref_start + 1];
    int snp_count[ref_end - ref_start + 1];
    int insert_count[ref_end - ref_start + 1];
    int delete_count[ref_end - ref_start + 1];
    vector< map<string, int> > AlleleFrequencyMap;
    vector< map<string, int> > AlleleFrequencyMapFwdStrand;
    vector< map<string, int> > AlleleFrequencyMapRevStrand;
    vector< set<string> > AlleleMap;

    // generate the image matrix of chunk_size (10kb) * feature_size (10)
    vector< vector<int> > image_matrix;

    image_matrix.resize(region_size + 1, vector<int>(feature_size));
    AlleleFrequencyMap.resize(region_size + 1);
    AlleleFrequencyMapFwdStrand.resize(region_size + 1);
    AlleleFrequencyMapRevStrand.resize(region_size + 1);
    AlleleMap.resize(region_size + 1);

    for (int i = 0; i < region_size + 1; i++) {
        for (int j = 0; j < feature_size; j++)
            image_matrix[i][j] = 0;
    }

    memset(coverage_vector, 0, sizeof(coverage_vector));
    memset(snp_count, 0, sizeof(snp_count));
    memset(insert_count, 0, sizeof(insert_count));
    memset(delete_count, 0, sizeof(delete_count));

    encode_reference_bases(image_matrix);

    // now iterate over all of the reads and populate the image matrix and coverage matrix
    for (auto &read:reads) {
        // this populates base_summaries and insert_summaries dictionaries
        if(print_colored_debug){
            cerr << BLUE << "Processing read: " << read.pos << " " <<read.query_name << RESET << endl;
        }
        if(read.mapping_quality > 0) {
            if (fahmid_check) {
                cerr << "Fahmid Check in gen sum: " << AlleleFrequencyMap.size() << endl;
            }

            populate_summary_matrix(image_matrix, coverage_vector, snp_count, insert_count, delete_count,
                                    AlleleFrequencyMap, AlleleFrequencyMapFwdStrand, AlleleFrequencyMapRevStrand, AlleleMap, read, min_snp_baseq, min_indel_baseq, train_mode);
        }
    }
  
    vector<long long> filtered_candidate_positions;
    bool snp_threshold_pass[ref_end - ref_start + 1];
    bool insert_threshold_pass[ref_end - ref_start + 1];
    bool delete_threshold_pass[ref_end - ref_start + 1];
    memset(snp_threshold_pass, 0, sizeof(snp_threshold_pass));
    memset(insert_threshold_pass, 0, sizeof(insert_threshold_pass));
    memset(delete_threshold_pass, 0, sizeof(delete_threshold_pass));

    // once the image matrix is generated, scale the counted values.
    for(int i=0;i<region_size;i++){
        double snp_fraction = snp_count[positions[i]-ref_start] / max(1.0, (double) coverage_vector[positions[i]-ref_start]);
        double insert_fraction = insert_count[positions[i]-ref_start] / max(1.0, (double) coverage_vector[positions[i]-ref_start]);
        double delete_fraction = delete_count[positions[i]-ref_start] / max(1.0, (double) coverage_vector[positions[i]-ref_start]);

        bool at_least_one_fraction_positive = (insert_fraction >0 || delete_fraction >0);
        bool indel_positive = (insert_fraction >= insert_freq_threshold || delete_fraction >= delete_freq_threshold);
        if(!indel_positive) {
            int scan_len = 20;
            // Instead of checking only one place for insert or delete, check a window of size 100bp centered around the candidate position
            if (i>scan_len && i<region_size-scan_len){
                for (int j=i-scan_len; j<i+scan_len; j++){
                    insert_fraction += (double) insert_count[positions[j]-ref_start] / max(1.0, (double) coverage_vector[positions[j]-ref_start]);
                    delete_fraction += (double) delete_count[positions[j]-ref_start] / max(1.0, (double) coverage_vector[positions[j]-ref_start]);
                }
            }
        }
        if(at_least_one_fraction_positive && print_colored_debug){
            cerr << GREEN << "At position: " << i << " positions[i]: " << positions[i] << " snp_fraction: " << snp_fraction << " insert_fraction: " << insert_fraction << " delete_fraction: " << delete_fraction 
            << " delete_count: "<< delete_count[positions[i]-ref_start] << " insert_count: " << insert_count[positions[i]-ref_start] <<" snp_count: " << snp_count[positions[i]-ref_start] << RESET << endl;
        }
        if(snp_fraction >= snp_freq_threshold || insert_fraction >= insert_freq_threshold || delete_fraction >= delete_freq_threshold) {
        
            // if(positions[i] >= candidate_region_start && positions[i] <= candidate_region_end && coverage_vector[positions[i]-ref_start] >= min_coverage_threshold) {
                filtered_candidate_positions.push_back(positions[i]);
                if(snp_fraction >= snp_freq_threshold) snp_threshold_pass[positions[i] - ref_start] = true;
                if(insert_fraction >= insert_freq_threshold) insert_threshold_pass[positions[i] - ref_start] = true;
                if(delete_fraction >= delete_freq_threshold) delete_threshold_pass[positions[i] - ref_start] = true;
            // }
        }

        for(int j=ImageOptionsRegion::BASE_INDEX_START; j < ImageOptionsRegion::BASE_INDEX_START + ImageOptionsRegion::BASE_INDEX_SIZE ; j++){
            if(image_matrix[i][j] >= 0)
                image_matrix[i][j] = (int) min(image_matrix[i][j], ImageOptionsRegion::MAX_COLOR_VALUE);
            else
                image_matrix[i][j] = (int) max(image_matrix[i][j], ImageOptionsRegion::MIN_COLOR_VALUE);
        }
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

    vector<CandidateImageSummary> all_candidate_images;
    // at this point all of the images are generated. So we can create the images for each candidate position.
    for(long long candidate_position : filtered_candidate_positions) {
        // if(print_colored_debug){
        //     cerr << BOLDCYAN << "Processing candidate position: " << candidate_position << " Relative pos: " << candidate_position - ref_start 
        //     << " AlleleMapSize: " << AlleleMap[candidate_position - ref_start].size() << RESET << endl;
        // }
        for (auto it=AlleleMap[candidate_position - ref_start].begin(); it!=AlleleMap[candidate_position - ref_start].end(); ++it) {
            CandidateImageSummary candidate_summary;
            candidate_summary.contig = contig;
            candidate_summary.position = candidate_position;
            bool debug = 0;
            if(debug) {
                cout << "-------------------------START----------------------------------------" << endl;
                cout << "Candidate position: " << candidate_position << endl;
                cout << "Coverage: " << coverage_vector[candidate_position - ref_start] << endl;
                // cout << "Candidates: " << endl;
            }

            candidate_summary.depth = min(coverage_vector[candidate_position-ref_start], ImageOptionsRegion::MAX_COLOR_VALUE);

            string candidate_string = *it;
            if(print_colored_debug && candidate_string[0]!='1'){ //skip snps for now
                cerr << BOLDCYAN << "Processing candidate position: " << candidate_position << RESET << endl;
                cerr << "Candidate string: " << candidate_string << endl;
            }
            bool is_deletion_end = false;
            if (candidate_string[0]=='E'){
                is_deletion_end = true;
            }
            int allele_depth = AlleleFrequencyMap[candidate_position - ref_start][candidate_string];
            int allele_depth_fwd = AlleleFrequencyMapFwdStrand[candidate_position - ref_start][candidate_string];
            int allele_depth_rev = AlleleFrequencyMapRevStrand[candidate_position - ref_start][candidate_string];
            double candidate_frequency = ((double) allele_depth / max(1.0, (double) candidate_summary.depth));
//            cout<<candidate_string<<" "<<allele_depth<<" "<<candidate_frequency<<endl;
            string candidate_allele = candidate_string.substr((is_deletion_end?2:1), candidate_string.length());
            // minimum 2 reads supporting the candidate or frequency is lower than 10
            
            if (allele_depth < candidate_support_threshold) {
                // cerr << "Skipping candidate at position: " << candidate_position << " as allele depth: " << allele_depth << " is less than candidate_support_threshold: " << candidate_support_threshold << endl;
                continue;
            }
            // see if candidate passes the candidate frequency threshold
            if (candidate_string[0] != '1' && candidate_frequency < indel_candidate_freq_threshold) {
                continue;
            }
            if (candidate_string[0] == '1' && candidate_frequency < snp_candidate_freq_threshold) {
                continue;
            }
            // If Candidate is SNP and we are skipping SNPs
            if ( candidate_string[0] == '1' and skip_snps == true) {
                continue;
            }
            // if Candidate is INDEL but we are skipping INDELs
            if ( candidate_string[0] != '1' and skip_indels == true) {
                continue;
            }
            // only pick type-specific candidates for each site
            if((candidate_string[0] == '1' && !snp_threshold_pass[candidate_position - ref_start]) ||
               (candidate_string[0] == '2' && !insert_threshold_pass[candidate_position - ref_start]) ||
               (candidate_string[(is_deletion_end?1:0)] == '3' && !delete_threshold_pass[candidate_position - ref_start])) {
                continue;
            }
            
            int base_index = (int) (candidate_position - ref_start + cumulative_observed_insert[candidate_position - ref_start]);
            char ref_base = reference_sequence[candidate_position - ref_start];
            if (train_mode) {
                vector<type_truth_record> hp1_truths = hp1_truth_alleles[base_index];
                vector<type_truth_record> hp2_truths = hp2_truth_alleles[base_index];
                int fixed_base_index = base_index;
                bool isIndel = candidate_string[0] != '1';
                if(hp1_truths.size()==0 && hp2_truths.size()==0 && isIndel) {
                    // This code is to handle the shift in the found alleles due to misplacement of the truth vcf
                    // It continually scans 150 bases on either side of the candidate position to find a truth allele
                    // Once found, it uses that position as the base_index for the candidate
                    // This base_index is later used to set the label for the candidate
                    fixed_base_index = -1;
                    for (int i = max(0, base_index - 300); i < min(base_index + 300, region_size); i++) {
                        int hap1_max_len = 0, hap2_max_len = 0;

                        if (hp1_truth_alleles[i].size() > 0 || hp2_truth_alleles[i].size() > 0){
                            if (candidate_string[0]=='2'){
                                for(const auto &rec: hp1_truth_alleles[i]){
                                    hap1_max_len = max(hap1_max_len, (int)rec.alt.length());
                                }
                                for (const auto &rec: hp2_truth_alleles[i]){
                                    hap2_max_len = max(hap2_max_len, (int)rec.alt.length());
                                }
                            }else{
                                for(const auto &rec: hp1_truth_alleles[i]){
                                    hap1_max_len = max(hap1_max_len, (int)rec.ref.length());
                                }
                                for (const auto &rec: hp2_truth_alleles[i]){
                                    hap2_max_len = max(hap2_max_len, (int)rec.ref.length());
                                }
                            
                            }
                        }
                        bool within_10pc_hap1 = candidate_string.length() >= 0.7 * hap1_max_len && candidate_string.length() <= 1.3 * hap1_max_len;
                        bool within_10pc_hap2 = candidate_string.length() >= 0.7 * hap2_max_len && candidate_string.length() <= 1.3 * hap2_max_len;
                        if (within_10pc_hap1 || within_10pc_hap2) {
                            fixed_base_index = i;
                            break;
                        }
                    }
                    if (fixed_base_index == -1) {
                        cerr << "No truth alleles found for candidate position: " << contig << ":" << candidate_position << endl;
                        fixed_base_index = base_index;
                        // continue;
                    }
                }
                hp1_truths = hp1_truth_alleles[fixed_base_index];
                hp2_truths = hp2_truth_alleles[fixed_base_index];
                vector<string> hp1_alleles;
                vector<string> hp2_alleles;
                for(const auto& truth_record : hp1_truths) {
                    if(truth_record.ref.length() > truth_record.alt.length()) {
                        //it's a delete
                        string alt_allele = truth_record.ref;
                        string ref_allele = truth_record.alt;
                        if(alt_allele.length() > 1 && ref_allele.length() > 1) {
                            int min_length = min(alt_allele.length(), ref_allele.length());
                            alt_allele = alt_allele.substr(0, alt_allele.length() - min_length + 1);
                        }
                        hp1_alleles.push_back(char(AlleleType::DELETE_ALLELE + '0') + alt_allele);
                    } else if(truth_record.ref.length() < truth_record.alt.length()) {
                        //it's an insert
                        string alt_allele = truth_record.alt;
                        string ref_allele = truth_record.ref;
                        if(alt_allele.length() > 1 && ref_allele.length() > 1) {
                            int min_length = min(alt_allele.length(), ref_allele.length());
                            alt_allele = alt_allele.substr(0, alt_allele.length() - min_length + 1);
                        }
                        hp1_alleles.push_back(char(AlleleType::INSERT_ALLELE + '0') + alt_allele);
                    } else if(truth_record.ref.length() == truth_record.alt.length()) {
                        //it's a SNP
                        string alt_allele = truth_record.alt;
                        string ref_allele = truth_record.ref;
                        if(alt_allele.length() > 1 && ref_allele.length() > 1) {
                            int min_length = min(alt_allele.length(), ref_allele.length());
                            alt_allele = alt_allele.substr(0, alt_allele.length() - min_length + 1);
                        }
                        hp1_alleles.push_back(char(AlleleType::SNP_ALLELE + '0') + alt_allele);
                    }
                }
                for(const auto& truth_record : hp2_truths) {
                    if(truth_record.ref.length() > truth_record.alt.length()) {
                        //it's a delete
                        string alt_allele = truth_record.ref;
                        string ref_allele = truth_record.alt;
                        if(alt_allele.length() > 1 && ref_allele.length() > 1) {
//                            cout<<"BEFORE: "<<alt_allele<<endl;
                            int min_length = min(alt_allele.length(), ref_allele.length());
                            alt_allele = alt_allele.substr(0, alt_allele.length() - min_length + 1);
//                            cout<<"AFTER: "<<alt_allele<<endl;
                        }

                        hp2_alleles.push_back(char(AlleleType::DELETE_ALLELE + '0') + alt_allele);
                    } else if(truth_record.ref.length() < truth_record.alt.length()) {
                        //it's an insert
                        string alt_allele = truth_record.alt;
                        string ref_allele = truth_record.ref;
                        if(alt_allele.length() > 1 && ref_allele.length() > 1) {
                            int min_length = min(alt_allele.length(), ref_allele.length());
                            alt_allele = alt_allele.substr(0, alt_allele.length() - min_length + 1);
                        }
                        hp2_alleles.push_back(char(AlleleType::INSERT_ALLELE + '0') + alt_allele);
                    } else if(truth_record.ref.length() == truth_record.alt.length()) {
                        //it's a SNP
                        string alt_allele = truth_record.alt;
                        string ref_allele = truth_record.ref;
                        if(alt_allele.length() > 1 && ref_allele.length() > 1) {
                            int min_length = min(alt_allele.length(), ref_allele.length());
                            alt_allele = alt_allele.substr(0, alt_allele.length() - min_length + 1);
                        }
                        hp2_alleles.push_back(char(AlleleType::SNP_ALLELE + '0') + alt_allele);
                    }
                }
                bool found_in_hp1 = false;
                bool found_in_hp2 = false;
//                cout<<"##############"<<endl;
//                cout<<"Candidate: "<<candidate_allele<<endl;
//                cout<<"HP1 truths:"<<endl;
                string cand_str = is_deletion_end ? candidate_string.substr(1, candidate_string.length()) : candidate_string;
                for(const auto& alt_hp1 : hp1_alleles) {
//                    cout<<alt_hp1<<endl;
                    if(alt_hp1 == candidate_string) {
                        found_in_hp1 = true;
                    }
                    // check if length of alt_hp1 is within 10% margin of candidate_allele
                    else if (isIndel && alt_hp1.length() >= 0.9 * cand_str.length() && alt_hp1.length() <= 1.1 * cand_str.length()) {
                        found_in_hp1 = true;
                    }
                    
                }

//                cout<<"HP2 truths:"<<endl;
                for(const auto& alt_hp2 : hp2_alleles) {
                //    cerr<<alt_hp2<<endl;
                    if(alt_hp2 == cand_str) {
                        found_in_hp2 = true;
                    }
                    else if (isIndel && alt_hp2.length() >= 0.9 * cand_str.length() && alt_hp2.length() <= 1.1 * cand_str.length()) {
                        found_in_hp2 = true;
                    }
                }
//                cout<<"##############"<<endl;

                int gt_label = 0;
                if(found_in_hp1 && found_in_hp2) {
                    gt_label = 2;
                } else if(found_in_hp1 || found_in_hp2){
                    gt_label = 1;
                }

                candidate_summary.base_label = labels[fixed_base_index];
                candidate_summary.type_label = gt_label;
                
                if(debug) {
                    cout << "BASE LABEL: " <<int(candidate_summary.base_label)<<endl;
                    cout << "TYPE LABEL: " <<int(candidate_summary.type_label)<<endl;
                    
                }
            } else {
                candidate_summary.base_label = 0;
                candidate_summary.type_label = 0;
            }
            
            int base_left = base_index - candidate_window_size / 2;
            int base_right = base_index + candidate_window_size / 2;
            // if(sifat_debug)
            //     cout<<"base left: "<<base_left<<" base right: "<<base_right<<" base_index: "<<base_index<<endl;
            // cout<<"--------------------------------------------------\n";
            candidate_summary.image_matrix.resize(candidate_window_size + 1, vector<int>(feature_size));
            // cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\t"<<endl;

            // cout << "FEATURE SIZE = " << image_matrix.size() << "----------------------------------------------------------\n\n";

            // now copy the entire feature matrix
            if (fahmid_check) {
                // cerr << "Fahmid Check" << endl;
            }

            for (int i = base_left; i <= base_right; i++) {
                for (int j = 0; j < feature_size; j++) {
                    if (i < 0 || i > region_size) {
                        candidate_summary.image_matrix[i - base_left][j] = 0;
                    } else {
                        candidate_summary.image_matrix[i - base_left][j] = image_matrix[i][j];
                    }
                }
            }
                        
            if (fahmid_check) {
                // cerr << "Fahmid Check" << endl;
            }

            // set type-specific feature values
            if(debug) {
                cout <<  "Candidate string:" <<candidate_string << ", allele_depth: " << allele_depth << ", candidate_summary depth: " << candidate_summary.depth << ", candidate_freq: " << candidate_frequency << endl;
                cout << "snp_thresh_pass: " <<snp_threshold_pass[candidate_position - ref_start] << ", ins_thresh_pass: " << insert_threshold_pass[candidate_position - ref_start] << ", del_thresh_pass: " << delete_threshold_pass[candidate_position - ref_start] << endl;
            }
            if (snp_threshold_pass[candidate_position - ref_start] && candidate_string[0] == '1') {
                if(debug) {
                    cout << "SNP " <<candidate_string << ",";
                    cout << AlleleFrequencyMap[candidate_position - ref_start][candidate_string] << endl;
                }
                // if(generate_negatives){
                //     allele_depth = 0;
                //     allele_depth_fwd = 0;
                //     allele_depth_rev = 0;
                // }
                int mid_index = candidate_window_size / 2;
                int forward_feature_index = get_feature_index(ref_base, candidate_string[1], false);
                int reverse_feature_index = get_feature_index(ref_base, candidate_string[1], true);
                candidate_summary.image_matrix[mid_index][1] = get_reference_feature_value(candidate_string[1]); // min((int) candidate_string.length() - 1, ImageOptionsRegion::MAX_COLOR_VALUE);
                candidate_summary.image_matrix[mid_index][5] = min(allele_depth_fwd, ImageOptionsRegion::MAX_COLOR_VALUE);
                candidate_summary.image_matrix[mid_index][16] = min(allele_depth_rev, ImageOptionsRegion::MAX_COLOR_VALUE);
                candidate_summary.image_matrix[mid_index][forward_feature_index] = (-1) * candidate_summary.image_matrix[mid_index][forward_feature_index];
                candidate_summary.image_matrix[mid_index][reverse_feature_index] = (-1) * candidate_summary.image_matrix[mid_index][reverse_feature_index];
                candidate_summary.candidates.push_back(candidate_string);
                candidate_summary.candidate_frequency.push_back(min(allele_depth, ImageOptionsRegion::MAX_COLOR_VALUE));
            } else if (insert_threshold_pass[candidate_position - ref_start] && candidate_string[0] == '2') {
                if(debug) {
                    cout << "INSERT" << " " << candidate_string << ",";
                    cout << AlleleFrequencyMap[candidate_position - ref_start][candidate_string] << endl;
                }
                int mid_index = candidate_window_size / 2;
                int forward_feature_index = get_feature_index(ref_base, 'I', false); // 12
                int reverse_feature_index = get_feature_index(ref_base, 'I', true); // 23
                candidate_summary.image_matrix[mid_index][2] = min((int) candidate_string.length() - 1, ImageOptionsRegion::MAX_COLOR_VALUE);
                candidate_summary.image_matrix[mid_index][6] = min(allele_depth_fwd, ImageOptionsRegion::MAX_COLOR_VALUE);
                candidate_summary.image_matrix[mid_index][17] = min(allele_depth_rev, ImageOptionsRegion::MAX_COLOR_VALUE);
                candidate_summary.image_matrix[mid_index][forward_feature_index] = (-1) * candidate_summary.image_matrix[mid_index][forward_feature_index];
                candidate_summary.image_matrix[mid_index][reverse_feature_index] = (-1) * candidate_summary.image_matrix[mid_index][reverse_feature_index];
                candidate_summary.candidates.push_back(candidate_string);
                candidate_summary.candidate_frequency.push_back(min(allele_depth, ImageOptionsRegion::MAX_COLOR_VALUE));
            } else if (delete_threshold_pass[candidate_position - ref_start] && candidate_string[(is_deletion_end?1:0)] == '3') {
                if(debug) {
                    cout << "DELETE: " << candidate_string << ",";
                    cout << AlleleFrequencyMap[candidate_position - ref_start][candidate_string] << endl;
                }
                int mid_index = candidate_window_size / 2;
                int del_len = (int) candidate_string.length() - 1 - (is_deletion_end?1:0);
                int end_index = min(mid_index + del_len - 1, candidate_window_size - 1);
                int forward_feature_index = get_feature_index(ref_base, 'D', false);
                int reverse_feature_index = get_feature_index(ref_base, 'D', true);
                candidate_summary.image_matrix[mid_index][3] = min((int) del_len, ImageOptionsRegion::MAX_COLOR_VALUE);
                candidate_summary.image_matrix[mid_index][7] = min(allele_depth_fwd, ImageOptionsRegion::MAX_COLOR_VALUE);
                candidate_summary.image_matrix[mid_index][18] = min(allele_depth_rev, ImageOptionsRegion::MAX_COLOR_VALUE);
                candidate_summary.image_matrix[mid_index][forward_feature_index] = (-1) * candidate_summary.image_matrix[mid_index][forward_feature_index];
                candidate_summary.image_matrix[mid_index][reverse_feature_index] = (-1) * candidate_summary.image_matrix[mid_index][reverse_feature_index];
                candidate_summary.candidates.push_back(candidate_string);
                candidate_summary.candidate_frequency.push_back(min(allele_depth, ImageOptionsRegion::MAX_COLOR_VALUE));
                // Update the image matrix for the rest of the window
                forward_feature_index = get_feature_index(ref_base, '*', false);
                reverse_feature_index = get_feature_index(ref_base, '*', true);
                if(!is_deletion_end){
                    for(int idx = mid_index + 1; idx <= end_index; idx++) {
                        candidate_summary.image_matrix[idx][3] = min((int) candidate_string.length() - 1, ImageOptionsRegion::MAX_COLOR_VALUE);
                        candidate_summary.image_matrix[idx][7] = min(allele_depth_fwd, ImageOptionsRegion::MAX_COLOR_VALUE);
                        candidate_summary.image_matrix[idx][18] = min(allele_depth_rev, ImageOptionsRegion::MAX_COLOR_VALUE);
                        candidate_summary.image_matrix[idx][forward_feature_index] = (-1) * candidate_summary.image_matrix[idx][forward_feature_index];
                        candidate_summary.image_matrix[idx][reverse_feature_index] = (-1) * candidate_summary.image_matrix[idx][reverse_feature_index];
                    }
                }else {
                    for(int idx = 0; idx < mid_index; idx++) {
                        candidate_summary.image_matrix[idx][3] = min((int) candidate_string.length() - 1, ImageOptionsRegion::MAX_COLOR_VALUE);
                        candidate_summary.image_matrix[idx][7] = min(allele_depth_fwd, ImageOptionsRegion::MAX_COLOR_VALUE);
                        candidate_summary.image_matrix[idx][18] = min(allele_depth_rev, ImageOptionsRegion::MAX_COLOR_VALUE);
                        candidate_summary.image_matrix[idx][forward_feature_index] = (-1) * candidate_summary.image_matrix[idx][forward_feature_index];
                        candidate_summary.image_matrix[idx][reverse_feature_index] = (-1) * candidate_summary.image_matrix[idx][reverse_feature_index];
                    }
                }
            }
            all_candidate_images.push_back(candidate_summary);
            if(debug) {
                debug_candidate_summary(candidate_summary, candidate_window_size, train_mode);
                cout << "-------------------------END----------------------------------------" << endl;
            }
        }
    }


    return all_candidate_images;
}


void RegionalSummaryGenerator::debug_print_matrix(vector<vector<int> > image_matrix, bool train_mode) {
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
    // for(int i=0; i < positions.size(); i++ ) {
    //     printf("%3lld\t", positions[i] % 100);
    // }
    cout << endl;
    int image_size = ImageOptionsRegion::REFERENCE_INDEX_SIZE + ImageOptionsRegion::SUPPORT_INDEX_SIZE + ImageOptionsRegion::BASE_INDEX_SIZE;
    for (int i = 0; i < image_size; i++) {
        cout<< ImageOptionsRegion::column_values[i] <<"\t";
        int region_size = (int) (ref_end - ref_start + total_observered_insert_bases + 1);

        for (int j = 0; j < region_size; j++) {
            printf("%3d\t", image_matrix[j][i]);
        }
        cout << endl;
    }

}


void RegionalSummaryGenerator::debug_candidate_summary(CandidateImageSummary candidate, int small_chunk_size, bool train_mode) {
    vector<string> decoded_base_lables {"RR", "RA", "RC", "RT", "RG", "R*", "R#", "AA", "AC", "AT", "AG", "A*", "A#", "CC", "CT", "CG", "C*", "C#", "TT", "TG", "T*", "T#", "GG", "G*", "G#", "**", "*#", "##"};
    vector<string> decoded_type_lables {"RR", "RS", "RI", "RD", "SS", "SI", "SD", "II", "ID", "DD" };
    cout << "------------- CANDIDATE PILEUP" << endl;
    cout<<"Contig: "<<candidate.contig<<endl;
    cout<<"Position: "<<candidate.position<<endl;
    cout<<"Type label: "<<(int)candidate.type_label<<" "<<decoded_type_lables[candidate.type_label]<<endl;
    cout<<"Base label: :"<<(int)candidate.base_label<<" "<<decoded_base_lables[candidate.base_label]<<endl;
    cout << "Refernce sequence: "<< reference_sequence << endl;

    long long candidate_ref_start = candidate.position - small_chunk_size / 2;
    long long candidate_ref_end = candidate.position + small_chunk_size / 2;
    cout << setprecision(3);
    for (long long i = candidate_ref_start; i <= candidate_ref_end; i++) {
        if (i == candidate_ref_start) cout << "POS:\t";
        // printf("%3lld\t", (i - candidate_ref_start) % 100);
        printf("%3lld\t", (i - candidate_ref_start));
    }
    cout << endl;

    for (long long i = candidate_ref_start; i <= candidate_ref_end; i++) {
        if(i==candidate_ref_start) cout<<"REF:\t";
        cout << "  " << reference_sequence[i - ref_start] << "\t";
        if (max_observed_insert[i - ref_start] > 0) {
            cout << max_observed_insert[i - ref_start] << "*\t";
            // for (uint64_t ii = 0; ii < max_observed_insert[i - ref_start]; ii++)cout << "  *" << "\t";
        }
    }
    cout << endl;


    if(train_mode) {
        for (long long i = candidate_ref_start; i <= candidate_ref_end; i++) {
            if(i==candidate_ref_start) cout<<"TRH1:\t";
            cout << "  " << labels_hp1[i - ref_start] << "\t";
            if (max_observed_insert[i - ref_start] > 0) {
                cout << max_observed_insert[i - ref_start] << "*\t";
                // for (uint64_t ii = 0; ii < max_observed_insert[i - ref_start]; ii++)cout << "  *" << "\t";
            }
        }
        cout << endl;

        for (long long i = candidate_ref_start; i <= candidate_ref_end; i++) {
            if(i==candidate_ref_start) cout<<"TRH2:\t";
            cout << "  " << labels_hp2[i - ref_start] << "\t";
            if (max_observed_insert[i - ref_start] > 0) {
                cout << max_observed_insert[i - ref_start] << "*\t";
                // for (uint64_t ii = 0; ii < max_observed_insert[i - ref_start]; ii++)cout << "  *" << "\t";
            }
        }
        cout << endl;

        for (long long i = candidate_ref_start; i <= candidate_ref_end; i++) {
            if(i==candidate_ref_start) cout<<"TRT1:\t";
            cout << "  " << variant_type_labels_hp1[i - ref_start] << "\t";
            if (max_observed_insert[i - ref_start] > 0) {
                cout << max_observed_insert[i - ref_start] << "*\t";
                // for (uint64_t ii = 0; ii < max_observed_insert[i - ref_start]; ii++)cout << "  *" << "\t";
            }
        }
        cout << endl;

        for (long long i = candidate_ref_start; i <= candidate_ref_end; i++) {
            if(i==candidate_ref_start) cout<<"TRT2:\t";
            cout << "  " << variant_type_labels_hp2[i - ref_start] << "\t";
            if (max_observed_insert[i - ref_start] > 0) {
                cout << max_observed_insert[i - ref_start] << "*\t";
                // for (uint64_t ii = 0; ii < max_observed_insert[i - ref_start]; ii++)cout << "  *" << "\t";
            }
        }
        cout << endl;
    }

    int image_size = candidate.image_matrix[0].size();
    for (int i = 0; i < image_size; i++) {
        cout<< ImageOptionsRegion::column_values[i] <<"\t";
        int region_size = candidate.image_matrix.size();

        for (int j = 0; j < region_size; j++) {
            printf("%3d\t", candidate.image_matrix[j][i]);
        }
        cout << endl;
    }

}