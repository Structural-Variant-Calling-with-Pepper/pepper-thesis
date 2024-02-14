import pysam
from pysam import VariantFile, VariantHeader
from pepper_variant.build import PEPPER_VARIANT
import collections
import math
import numpy as np


class VCFWriter:
    def __init__(self, reference_file_path, sample_name, output_dir, filename):
        self.fasta_handler = PEPPER_VARIANT.FASTA_handler(reference_file_path)
        contigs = self.fasta_handler.get_chromosome_names()
        self.contigs = contigs
        self.vcf_header = self.get_vcf_header(sample_name, contigs)
        self.output_dir = output_dir
        self.vcf_file_name = self.output_dir + filename + '.vcf.gz'
        self.vcf_file = VariantFile(self.vcf_file_name, 'w', header=self.vcf_header)

    def __del__(self):
        # close the files
        self.vcf_file.close()

        # index the files
        pysam.tabix_index(self.vcf_file_name, preset="vcf", force=True)

    def dump_candidates(self, candidates):
        for candidate in candidates:
            contig, ref_start, ref_end, ref_allele, alts, genotype, depth, variant_allele_support = candidate

            ## This is where things become tricky
            ## Change or comment this block and truvari doesn't work anymore
            # Example region: 20:286892-286931
            if len(alts) > 1:
                alts = [alts[0]]
                variant_allele_support = [variant_allele_support[0]]

            qual = 10
            alleles = tuple([ref_allele]) + tuple(alts)
            vcf_record = self.vcf_file.new_record(contig=str(contig), start=ref_start,
                                                  stop=ref_end, id='.', qual=qual,
                                                  filter='PASS', alleles=alleles, GT=genotype,
                                                  GQ=qual, DP=depth, AD=variant_allele_support)

            self.vcf_file.write(vcf_record)

    def get_vcf_header(self, sample_name, contigs):
        header = VariantHeader()

        items = [('ID', "PASS"),
                 ('Description', "All filters passed")]
        header.add_meta(key='FILTER', items=items)
        items = [('ID', "refCall"),
                 ('Description', "Call is homozygous")]
        header.add_meta(key='FILTER', items=items)
        items = [('ID', "lowGQ"),
                 ('Description', "Low genotype quality")]
        header.add_meta(key='FILTER', items=items)
        items = [('ID', "lowQUAL"),
                 ('Description', "Low variant call quality")]
        header.add_meta(key='FILTER', items=items)
        items = [('ID', "conflictPos"),
                 ('Description', "Overlapping record")]
        header.add_meta(key='FILTER', items=items)
        items = [('ID', "GT"),
                 ('Number', 1),
                 ('Type', 'String'),
                 ('Description', "Genotype")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "DP"),
                 ('Number', 1),
                 ('Type', 'Integer'),
                 ('Description', "Depth")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "AD"),
                 ('Number', "A"),
                 ('Type', 'Integer'),
                 ('Description', "Allele depth")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "VAF"),
                 ('Number', "A"),
                 ('Type', 'Float'),
                 ('Description', "Variant allele fractions.")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "AP"),
                 ('Number', "A"),
                 ('Type', 'Float'),
                 ('Description', "Maximum variant allele probability for each allele.")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "GT"),
                 ('Number', 1),
                 ('Type', 'String'),
                 ('Description', "Genotype")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "GQ"),
                 ('Number', 1),
                 ('Type', 'Float'),
                 ('Description', "Genotype Quality")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "REP"),
                 ('Number', 1),
                 ('Type', 'String'),
                 ('Description', "If set to 1 then variant site is considered to be ina LowCompexity repeat region")]
        header.add_meta(key='FORMAT', items=items)

        sqs = self.fasta_handler.get_chromosome_names()
        for sq in sqs:
            if sq not in contigs:
                continue
            sq_id = sq
            ln = self.fasta_handler.get_chromosome_sequence_length(sq)
            header.contigs.add(sq_id, length=ln)

        header.add_sample(sample_name)

        return header
