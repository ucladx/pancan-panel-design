#!/usr/bin/env python3
'''
@description : Given a gnomAD VCF, shortlist SNPs for a capture-based panel to detect allelic imbalance
@created : 09/16/2020
@author : Cyriac Kandoth

'''
from __future__ import division, print_function
import os, sys, argparse, time, pysam

def main():

    parser = argparse.ArgumentParser(prog='select_gnomad_snps.py', description='shortlist GRCh38 SNPs from gnomAD v3, best at detecting allelic imbalance', usage='python3 %(prog)s [options]')
    parser.add_argument("--gnomad-vcf", action="store", metavar="VCF", dest="gnomad_vcf", required=True, type=str, help="gnomAD v3 VCF, bgzipped and indexed")
    parser.add_argument("--max-snps", action="store", metavar="NUM", dest="max_snps", required=False, type=int, default=50000, help="no more than this many SNPs will be shortlisted")
    parser.add_argument("--output-bed", action="store", metavar="BED", dest="output_bed", required=True, type=str, help="path to save resulting SNP BED file")
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.gnomad_vcf)
    shortlist = dict()
    for var in vcf.fetch():
        # gnomAD v3 lists multiallelic positions in separate rows, so we can safely select only the first ALT
        alt = var.alts[0]
        # skip rows not containing snps, FILTER!=PASS, MQ<50, or are depleted for homozygotes
        if(len(var.ref) > 1 or len(alt) > 1 or 'PASS' not in var.filter.keys() or var.info.get('MQ') < 50 or var.info.get('nhomalt')[0] == 0):
            continue
        # count how many of 9 subpopulations have a frequency of heterozygotes >25%
        min_het_freq = 0.25
        het_subpops = 0
        sum_het_freq = 0
        for subpop in ["ami", "oth", "afr", "sas", "asj", "fin", "amr", "nfe", "eas"]:
            alt_count = var.info.get('AC_' + subpop)[0]
            hom_alt_count = var.info.get('nhomalt_' + subpop)[0]
            total_alleles = var.info.get('AN_' + subpop)
            het_freq = 0
            if(alt_count > 0 and total_alleles > 0):
                het_freq = (alt_count - 2*hom_alt_count) / (total_alleles/2)
            if(het_freq > min_het_freq):
                het_subpops += 1
            sum_het_freq += het_freq
        avg_het_freq = sum_het_freq / 9
        var_bed = '\t'.join(str(x) for x in [var.chrom, var.pos-1, var.pos, '.' if var.id is None else var.id])
        # select snps meeting the minimum frequency of heterozygotes in all 9 subpopulations
        if(het_subpops >= 9 and var_bed not in shortlist):
            shortlist[var_bed] = avg_het_freq

    with open(args.output_bed, 'w') as bed:
        snp_count = 0
        for var_bed, avg_het_freq in sorted(shortlist.items(), key=lambda x: x[1], reverse=True):
            print('\t'.join([var_bed, str(round(avg_het_freq, 6)), "+"]), file=bed)
            snp_count += 1
            if(snp_count >= args.max_snps):
                break

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    total_time = end_time - start_time
    sys.stderr.write("Runtime: %0.2f seconds\n" % total_time)
