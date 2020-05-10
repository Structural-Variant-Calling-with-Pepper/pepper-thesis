import argparse
from collections import defaultdict
from pysam import VariantFile
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def plot_distributions(vafs, true_vafs, false_vafs):
    sns.set(rc={"figure.figsize": (20, 10)})
    sns.set_style("white")
    # plt.hist(vafs, bins=100, color='blue', alpha=0.4)
    plt.hist([true_vafs, false_vafs], bins=100, density=False, histtype='bar', color=['green', 'red'], alpha=0.4, stacked=True, label=['True variants', 'False positives'])
    # plt.hist(false_vafs, bins=100, color='red', alpha=0.4, stacked=True)
    plt.xlim((0.00, 1.10))
    # plt.ylim((0, 500))
    plt.legend(fontsize='x-large')
    # plt.xticks(np.arange(0, 1, step=0.10), fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel("Allele frequency", fontsize='24')
    plt.ylabel("Count", fontsize='24')

    plt.title("Stacked histogram showing TP and FP distribution at different frequency intervals.", fontsize='20')
    # plt.show()
    # exit()
    output_file_name = "./VAF_distribution.png"
    plt.savefig(output_file_name, format='png', dpi=300, quality=95)


def plot_distributions_q(true_qvals, false_qvals):
    sns.set(rc={"figure.figsize": (20, 10)})
    sns.set_style("white")
    # plt.hist(vafs, bins=100, color='blue', alpha=0.4)
    plt.hist([true_qvals, false_qvals], bins=100, density=False, histtype='bar', color=['green', 'red'], alpha=0.4, stacked=True, label=['True variants', 'False positives'])
    # plt.hist(false_vafs, bins=100, color='red', alpha=0.4, stacked=True)
    plt.xlim((0.00, 1.10))
    # plt.ylim((0, 500))
    plt.legend(fontsize='x-large')
    # plt.xticks(np.arange(0, 1, step=0.10), fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel("P-value", fontsize='24')
    plt.ylabel("Count", fontsize='24')

    plt.title("Stacked histogram showing TP and FP distribution at different P-value intervals.", fontsize='20')
    # plt.show()
    # exit()
    output_file_name = "./P_val_distribution.png"
    plt.savefig(output_file_name, format='png', dpi=300, quality=95)


def plot_vaf_and_q_2d(true_allele_q_vaf, false_allele_q_vaf, q_vaf):
    sns.set(rc={"figure.figsize": (20, 10)})
    sns.set_style("white")
    q, vaf, pred = zip(*q_vaf)
    plt.scatter(vaf, q, c=pred)
    plt.xlabel("VAF", fontsize='24')
    plt.ylabel("Q-value", fontsize='24')
    plt.xlim((0.00, 0.25))
    plt.show()


def calculate_stats(truth_vcf, vcf):
    untagged_vcf = VariantFile(vcf)
    positional_norm_q = defaultdict()
    for rec in untagged_vcf.fetch():
        positional_norm_q[rec.pos] = rec.qual

    vcf_in1 = VariantFile(truth_vcf)

    vafs_of_true_alleles = list()
    vafs_of_false_alleles = list()
    q_of_true_alleles = list()
    q_of_false_alleles = list()
    total_alts = 0
    total_true_calls = 0
    total_false_calls = 0
    true_allele_q_vaf = list()
    false_allele_q_vaf = list()
    q_vaf = list()

    VAF_Threshold = 1.10

    all_allele_frequencies = list()
    total_recs = 0
    total_multi_allelic_sites = 0
    for rec in vcf_in1.fetch():
        total_recs += 1
        # ['__class__', '__delattr__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__setstate__', '__sizeof__', '__str__', '__subclasshook__',
        # 'alleles', 'alts', 'chrom', 'contig', 'copy', 'filter', 'format', 'header', 'id', 'info', 'pos', 'qual', 'ref', 'rid', 'rlen', 'samples', 'start', 'stop', 'translate']
        alts = rec.alts
        total_alts += len(alts)
        for sample in rec.samples:
            if 'VAF' in rec.samples[sample].keys():
                vafs = list(rec.samples[sample]['VAF'])
            else:
                vafs = [0.0] * len(alts)

            gts = list(rec.samples[sample]['GT'])

            true_index = []
            for gt in gts:
                if gt != 0:
                    true_index.append(gt-1)

            if true_index:
                # if positional_norm_q[rec.pos] == 0:
                #     print(rec)
                if min(vafs) <= VAF_Threshold:
                    q_of_true_alleles.append(positional_norm_q[rec.pos])
            else:
                if min(vafs) <= VAF_Threshold:
                    q_of_false_alleles.append(positional_norm_q[rec.pos])

            if len(alts) > 1:
                total_multi_allelic_sites += 1
            for i, (alt, vaf) in enumerate(zip(alts, vafs)):
                if i in true_index:
                    vafs_of_true_alleles.append(vaf)
                    true_allele_q_vaf.append((positional_norm_q[rec.pos], min(1.1, vaf)))

                    if vaf <= VAF_Threshold:
                        q_vaf.append((positional_norm_q[rec.pos], min(1.1, vaf), 'Green'))

                    total_true_calls += 1
                else:
                    vafs_of_false_alleles.append(vaf)
                    false_allele_q_vaf.append((positional_norm_q[rec.pos], min(1.1, vaf)))

                    if vaf <= VAF_Threshold:
                        q_vaf.append((positional_norm_q[rec.pos], min(1.1, vaf), 'Red'))

                    total_false_calls += 1

            for vaf in vafs:
                all_allele_frequencies.append(round(vaf, 3))

    print("Records:\t", total_recs)
    print("Multi-alleleic:\t", total_false_calls, "(" + str(int(100 * (total_false_calls/total_alts))) + "%)")
    print("Alt alleles:\t", total_alts)
    print("True alleles:\t", total_true_calls, "(" + str(int(100 * (total_true_calls/total_alts))) + "%)")
    print("False alleles:\t", total_false_calls, "(" + str(int(100 * (total_false_calls/total_alts))) + "%)")

    plot_distributions(all_allele_frequencies, vafs_of_true_alleles, vafs_of_false_alleles)
    # plot_distributions_q(q_of_true_alleles, q_of_false_alleles)
    # plot_vaf_and_q_2d(true_allele_q_vaf, false_allele_q_vaf, q_vaf)



def add_merge_vcf_arguments(parser):
    parser.add_argument(
        "-vt",
        "--truth_vcf",
        type=str,
        required=True,
        help="VCF of haplotype 1."
    )
    parser.add_argument(
        "-v",
        "--vcf",
        type=str,
        required=True,
        help="VCF of haplotype 1."
    )
    return parser


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="PEPPER is a RNN based polisher for polishing ONT-based assemblies. "
                                                 "It works in three steps:\n"
                                                 "1) make_images: This module takes alignment file and coverts them"
                                                 "to HDF5 files containing summary statistics.\n"
                                                 "2) run_inference: This module takes the summary images and a"
                                                 "trained neural network and generates predictions per base.\n"
                                                 "3) find_snps: This module takes the inference files as input and "
                                                 "finds possible SNP sites.\n",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--version",
        default=False,
        action='store_true',
        help="Show version."
    )
    add_merge_vcf_arguments(parser)
    FLAGS, unparsed = parser.parse_known_args()
    calculate_stats(FLAGS.truth_vcf, FLAGS.vcf)