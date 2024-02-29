from pepper_variant.build import PEPPER_VARIANT
import numpy as np
from pysam import VariantFile
from pepper_variant.modules.python.Options import ImageSizeOptions, AlingerOptions, ConsensCandidateFinder
import sys
from colorama import Fore, Style

fahmid_check = False

class AlignmentSummarizer:
    """

    """
    def __init__(self, bam_handler, fasta_handler, chromosome_name, region_start, region_end):
        if fahmid_check :
            sys.stderr.write(Fore.RED + f"Fahmid Check in Allignement Summarizer Constructor, region start: {region_start} and region end: {region_end}\n" + Style.RESET_ALL)

        self.bam_handler = bam_handler
        self.fasta_handler = fasta_handler
        self.chromosome_name = chromosome_name
        self.region_start_position = region_start
        self.region_end_position = region_end
        self.print_colored_debug = False

    @staticmethod
    def range_intersection_bed(interval, bed_intervals, options):
        left = interval[0]
        right = interval[1]

        intervals = []
        for i in range(0, len(bed_intervals)):
            bed_left = bed_intervals[i][0]
            bed_right = bed_intervals[i][1]

            if bed_right < left:
                continue
            elif bed_left > right:
                continue
            else:
                if not options.train_mode:
                    left_bed = max(left, bed_left)
                    right_bed = min(right, bed_right)
                else:
                    left_bed = bed_left
                    right_bed = bed_right
                    if i == 0 and left_bed < 300:
                        left_bed = 301
                    elif i == len(bed_intervals) - 1:
                        right_bed -= 301
                
                intervals.append([left_bed, right_bed])

        return intervals

    def get_truth_vcf_records(self, vcf_file, region_start, region_end):
        truth_vcf_file = VariantFile(vcf_file)
        all_records = truth_vcf_file.fetch(self.chromosome_name, region_start, region_end)
        haplotype_1_records = []
        haplotype_2_records = []
        for record in all_records:
            # filter only for PASS variants
            if 'PASS' not in record.filter.keys():
                continue

            genotype = []
            for sample_name, sample_items in record.samples.items():
                sample_items = sample_items.items()
                if self.print_colored_debug:
                    sys.stderr.write(Fore.GREEN + f"Sample name : {sample_name}\n" + Fore.RESET)
                    sys.stderr.write(Fore.GREEN + f"Sample items : {sample_items}\n" + Fore.RESET)
                for name, value in sample_items:
                    if name == 'GT':
                        genotype = value

            for hap, alt_location in enumerate(genotype):
                if self.print_colored_debug:
                    sys.stderr.write(Fore.GREEN + f"Hap : {hap}, Alt location : {alt_location}\n" + Fore.RESET)
                if alt_location == 0:
                    continue
                if hap == 0:
                    truth_variant = PEPPER_VARIANT.type_truth_record(record.contig, record.start, record.stop, record.alleles[0], record.alleles[alt_location])
                    haplotype_1_records.append(truth_variant)
                else:
                    truth_variant = PEPPER_VARIANT.type_truth_record(record.contig, record.start, record.stop, record.alleles[0], record.alleles[alt_location])
                    haplotype_2_records.append(truth_variant)

        return haplotype_1_records, haplotype_2_records
    
    def get_previous_interval(self, bed_intervals, interval):
        for i in range(0, len(bed_intervals)):
            bed_left = bed_intervals[i][0]
            bed_right = bed_intervals[i][1]

            if bed_left == interval[0] and bed_right == interval[1]:
                if i == 0:
                    return None
                else:
                    return bed_intervals[i-1][1]
        return None

    def create_summary(self, options, bed_list, thread_id):
        """
        Core process to generate image summaries.
        :param options: Options for image generation
        :param bed_list: List of regions from a bed file. [GIAB high-confidence region]
        :param thread_id: Process id.
        :return:
        """
        all_candidate_images = []
       
        truth_regions = []
        # now intersect with bed file
        if bed_list is not None:
            intersected_truth_regions = []
            if self.chromosome_name in bed_list.keys():
                reg = AlignmentSummarizer.range_intersection_bed([self.region_start_position, self.region_end_position], bed_list[self.chromosome_name], options)
                intersected_truth_regions.extend(reg)
            truth_regions = intersected_truth_regions
        if options.train_mode:
            if not truth_regions:
                # sys.stderr.write(TextColor.GREEN + "INFO: " + log_prefix + " NO TRAINING REGION FOUND.\n"
                #                  + TextColor.END)
                return None
            if self.print_colored_debug:
                sys.stderr.write(Fore.RED + f"Provided: {self.chromosome_name}:{self.region_start_position}-{self.region_end_position}\n" + Style.RESET_ALL)
                sys.stderr.write(Fore.RED + f"Truth deletion: {str(truth_regions)}\n" + Fore.RESET)
            # print("Higher Region : ( "+str(self.region_start_position)+" , "+str(self.region_end_position)+" )")
            i = 0
            # for item in bed_list[self.chromosome_name]:
            #     # check if item is within higher region
            #     if item[1] < self.region_start_position or item[0] > self.region_end_position:
            #         continue
            #     sys.stderr.write(f"item : {i}, {item}\n")
            #     i +=1
            for region in truth_regions:
                sub_region_start = region[0]-ConsensCandidateFinder.REGION_SAFE_BASES
                sub_region_end = region[1] + ConsensCandidateFinder.REGION_SAFE_BASES
                # sys.stderr.write("Sub-Higher Region : ( "+str(region[0])+" , "+str(region[1])+" )\n")

                # --------------------------------------for positiv sample ----------------------------------------------
                if self.print_colored_debug:
                    sys.stderr.write(Fore.YELLOW + f"Padding: {ConsensCandidateFinder.REGION_SAFE_BASES}\n" + Fore.RESET)
                    sys.stderr.write(Fore.YELLOW + f"Window: {sub_region_start}-{sub_region_end}\n" + Fore.RESET)
                # --------------------------------------for negative sample ---------------------------------------------
                # curr_interval_start = region[0]
                # prev_interval_end = self.get_previous_interval(bed_list[self.chromosome_name],region)

                # if prev_interval_end is None:
                #     continue
                # else:
                #     print("neg start : ",prev_interval_end,"  end : ",curr_interval_start) 

                # jump_size = 10000
                # pad  = 1000
                # for st in range(prev_interval_end+pad,curr_interval_start-pad-jump_size,jump_size):
                        

                    # window_start = st 
                    # window_end = st + jump_size

                # ===============================================================================================================


                region_start = max(0, sub_region_start)
                region_end = sub_region_end

                all_reads = self.bam_handler.get_reads(self.chromosome_name,
                                                        region_start,
                                                        region_end + 1,
                                                        options.include_supplementary,
                                                        options.min_mapq,
                                                        options.min_snp_baseq)

                
        
                total_reads = len(all_reads)
                total_allowed_reads = int(min(AlingerOptions.MAX_READS_IN_REGION, options.downsample_rate * total_reads))
                if self.print_colored_debug:
                    sys.stderr.write(Fore.YELLOW + f"Total reads: {total_reads}, Total allowed reads: {total_allowed_reads}\n" + Fore.RESET)
                # print("Total reads: ", total_reads)
                # print("Total allowed reads: ", total_allowed_reads)
                if total_reads > total_allowed_reads:
                    # https://github.com/google/nucleus/blob/master/nucleus/util/utils.py
                    # reservoir_sample method utilized here
                    random = np.random.RandomState(AlingerOptions.RANDOM_SEED)
                    sample = []
                    for i, read in enumerate(all_reads):
                        if len(sample) < total_allowed_reads:
                            sample.append(read)
                        else:
                            j = random.randint(0, i + 1)
                            if j < total_allowed_reads:
                                sample[j] = read
                    all_reads = sample

                total_reads = len(all_reads)

                if total_reads == 0:
                    continue
                # pstr = " start" if i==0 else " end"
                # print("\t\t\t\t\t\t\tSub Region Pos: "+str(region[i])+pstr)
                # print("\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tSub Region : ( "+str(region_start)+" , "+str(region_end)+" )")
                # get vcf records from truth
                truth_hap1_records, truth_hap2_records = self.get_truth_vcf_records(options.truth_vcf, region_start, region_end)
                if self.print_colored_debug:
                    sys.stderr.write(Fore.BLUE + f"Iteration: {i}\n" + Fore.RESET)
                    sys.stderr.write(Fore.GREEN + f"Truth hap1 records : {truth_hap1_records}\n" + Fore.RESET)
                    sys.stderr.write(Fore.GREEN + f"Truth hap2 records : {truth_hap2_records}\n" + Fore.RESET)
                # ref_seq should contain region_end_position base
                ref_seq = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                                    region_start,
                                                                    region_end + 1)
                # sys.stderr.write(Fore.GREEN + f"Ref seq : {ref_seq}\n" + Fore.RESET)
                regional_summary = PEPPER_VARIANT.RegionalSummaryGenerator(self.chromosome_name, region_start, region_end, ref_seq)

                regional_summary.generate_max_insert_summary(all_reads)

                regional_summary.generate_labels(truth_hap1_records, truth_hap2_records)
                
                if self.print_colored_debug:
                    sys.stderr.write(Fore.LIGHTWHITE_EX + f"Regional summary : {regional_summary}\n" + Fore.RESET)
                
                if fahmid_check :
                    sys.stderr.write(Fore.RED + f"Fahmid Check in Allignement Summarizer Just b4 gen sum, region start: {self.region_start_position} and region end: {self.region_end_position}\n" + Style.RESET_ALL)

                candidate_image_summary = regional_summary.generate_summary(all_reads,
                                                                            options.min_snp_baseq,
                                                                            options.min_indel_baseq,
                                                                            options.snp_frequency,
                                                                            options.insert_frequency,
                                                                            options.delete_frequency,
                                                                            options.min_coverage_threshold,
                                                                            options.snp_candidate_frequency_threshold,
                                                                            options.indel_candidate_frequency_threshold,
                                                                            options.candidate_support_threshold,
                                                                            options.skip_indels,
                                                                            self.region_start_position,
                                                                            self.region_end_position,
                                                                            ImageSizeOptions.CANDIDATE_WINDOW_SIZE,
                                                                            ImageSizeOptions.IMAGE_HEIGHT,
                                                                            options.train_mode)
                if self.print_colored_debug:
                    sys.stderr.write(Fore.LIGHTWHITE_EX + f"Candidate image summary : {candidate_image_summary}\n" + Fore.RESET)
                total_ref_examples = 0
                for candidate in candidate_image_summary:
                    if candidate.type_label == 0:
                        total_ref_examples += 1

                picked_refs = 0
                random_sampling = np.random.uniform(0.0, 1.0, total_ref_examples)
                random_sampling_index = 0
                for candidate in candidate_image_summary:
                    if candidate.type_label == 0:
                        random_draw = random_sampling[random_sampling_index]
                        random_sampling_index += 1
                        if random_draw <= options.random_draw_probability:
                            all_candidate_images.append(candidate)
                            picked_refs += 1
                    else:
                        all_candidate_images.append(candidate)
           
        else:
            ConsensCandidateFinder.REGION_SAFE_BASES = 0
            if truth_regions is None:
                truth_regions = [[self.region_start_position, self.region_end_position]]
            for region in truth_regions:
                region_start = max(0, region[0] - ConsensCandidateFinder.REGION_SAFE_BASES)
                region_end = region[1] + ConsensCandidateFinder.REGION_SAFE_BASES
                # sys.stderr.write(Fore.YELLOW + f"Region start: {region_start}, Region end: {region_end}\n" + Fore.RESET)
                all_reads = self.bam_handler.get_reads(self.chromosome_name,
                                                    region_start,
                                                    region_end,
                                                    options.include_supplementary,
                                                    options.min_mapq,
                                                    options.min_snp_baseq)

                total_reads = len(all_reads)
                total_allowed_reads = int(min(AlingerOptions.MAX_READS_IN_REGION, options.downsample_rate * total_reads))
                # print("Total reads: ", total_reads)
                # print("Total allowed reads: ", total_allowed_reads)
                if total_reads > total_allowed_reads:
                    # sys.stderr.write("INFO: " + log_prefix + "HIGH COVERAGE CHUNK: " + str(total_reads) + " Reads.\n")
                    # https://github.com/google/nucleus/blob/master/nucleus/util/utils.py
                    # reservoir_sample method utilized here
                    random = np.random.RandomState(AlingerOptions.RANDOM_SEED)
                    sample = []
                    for i, read in enumerate(all_reads):
                        if len(sample) < total_allowed_reads:
                            sample.append(read)
                        else:
                            j = random.randint(0, i + 1)
                            if j < total_allowed_reads:
                                sample[j] = read
                    all_reads = sample

                total_reads = len(all_reads)

                if total_reads == 0:
                    return None

                # ref_seq should contain region_end_position base
                ref_seq = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                                    region_start,
                                                                    region_end + 1)

                regional_summary = PEPPER_VARIANT.RegionalSummaryGenerator(self.chromosome_name, region_start, region_end, ref_seq)
                regional_summary.generate_max_insert_summary(all_reads)

                if fahmid_check :
                    sys.stderr.write(Fore.RED + f"Fahmid Check in Allignement Summarizer Just b4 gen sum, region start: {self.region_start_position} and region end: {self.region_end_position}\n" + Style.RESET_ALL)

                candidate_image_summary = regional_summary.generate_summary(all_reads,
                                                                            options.min_snp_baseq,
                                                                            options.min_indel_baseq,
                                                                            options.snp_frequency,
                                                                            options.insert_frequency,
                                                                            options.delete_frequency,
                                                                            options.min_coverage_threshold,
                                                                            options.snp_candidate_frequency_threshold,
                                                                            options.indel_candidate_frequency_threshold,
                                                                            options.candidate_support_threshold,
                                                                            options.skip_indels,
                                                                            self.region_start_position,
                                                                            self.region_end_position,
                                                                            ImageSizeOptions.CANDIDATE_WINDOW_SIZE,
                                                                            ImageSizeOptions.IMAGE_HEIGHT,
                                                                            options.train_mode)

                all_candidate_images.extend(candidate_image_summary)

        return all_candidate_images