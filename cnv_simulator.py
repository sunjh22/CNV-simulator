#! /usr/bin/env python

import argparse
import os
from genome_simu import *


class CapitalisedHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        if prefix is None:
            prefix = 'Usage: '
            return super(CapitalisedHelpFormatter, self).add_usage(usage, actions, groups, prefix)

def main():
    parser = argparse.ArgumentParser(add_help=True, formatter_class=CapitalisedHelpFormatter, \
                                     description="Simulate whole-genome level haplotype-resolved CNVs and generate NGS reads")
    parser._positionals.title = 'Positional arguments'
    parser._optionals.title = 'Optional arguments'
    parser.add_argument("-v", "--version", action="version", version="CNV-simulator v0.1", \
                        help="Show program's version number and exit.")

    parser.add_argument("genome", type=str, help="reference genome file in fasta format")
    parser.add_argument("access", type=str, help="genomic accessible region file")

    parser.add_argument("-o", "--output_dir", type=str, default="tmp_simulaiton", \
                        help="a name to be used to create the output directory (overrides existing directory with the same name).")
    parser.add_argument("-a" ,"--prefix", type=str, default="sample", \
                        help="prefix of simulated genome and fastq files")
    parser.add_argument("-l", "--read_length", type=int, default=150, \
                        help="read length (bp)")
    parser.add_argument("-i" ,"--cnv_list", type=str, default=None, \
                        help="path to a CNV list file in BED format chr | start | end | cn. If not passed, it is randomly generated using CNV list parameters below")
    parser.add_argument("-A", "--append_list", type=str, default='False', \
                        help="generate some new CNVs and append them into the user-provided CNV list")
    parser.add_argument("-c" ,"--coverage", type=float, default=1, \
                        help="coverage for reads simulation")
    
    cnv_sim_group = parser.add_argument_group('CNV list parameters', "parameters to be used to generate CNV list if CNV list is not passed")

    cnv_sim_group.add_argument("-n", "--cnv_number", type=int, default=200, \
                               help="the number of CNVs to be simulated")
    cnv_sim_group.add_argument("-b", "--cnv_minimum_length", type=int, default=5000, \
                               help="minimum length of each CNV region")
    cnv_sim_group.add_argument("-B", "--cnv_maximum_length", type=int, default=2000000, \
                               help="maximum length of each CNV region")
    cnv_sim_group.add_argument("-e", "--cnv_expo_scale", type=int, default=200000, \
                               help="scale of a exponential distribution for simulating the length of CNV")
    cnv_sim_group.add_argument("-p", "--amplify_prop", type=float, default=0.50, \
                        help="percentage of amplifications in range [0.0: 1.0].")
    
    args = parser.parse_args()

    genomeFile = args.genome
    accessFile = args.access
    outDir = args.output_dir
    prefix = args.prefix
    readLength = args.read_length
    coverage = args.coverage / 2
    cnvNum = args.cnv_number
    minCnvLength = args.cnv_minimum_length
    maxCnvLength = args.cnv_maximum_length
    expCnvLength = args.cnv_expo_scale
    amplifyProp = args.amplify_prop

    paternalGenomeFile = os.path.join(outDir, prefix+'_paternal.fa')
    maternalGenomeFile = os.path.join(outDir, prefix+'_maternal.fa')
    paternalReadsFile = os.path.join(outDir, prefix+'_paternal_')
    maternalReadsFile = os.path.join(outDir, prefix+'_maternal_')

    if not os.path.exists(outDir):
        os.makedirs(outDir, exist_ok=True)

    if args.cnv_list is not None:
        cnvList = args.cnv_list
        log("An input CNV list is detected at " + repr(cnvList))
        if args.append_list == 'True':
            log("Appending mode is activated, simulate new CNVs and append them into existing list")
            GenerateCNVList(access_file=accessFile, output_file=cnvList, cnv_num=cnvNum, p_amplify=amplifyProp, min_length=minCnvLength, max_length=maxCnvLength, exp_length=expCnvLength, append=True)
            subprocess.call(['sort', '-Vk 1 -k 2,3n', cnvList, '-o', cnvList])
    else:
        log("No CNV list was detected, start simulating a CNV list")
        cnvList = os.path.join(outDir, prefix+'_cnvList.bed')
        GenerateCNVList(access_file=accessFile, output_file=cnvList, cnv_num=cnvNum, p_amplify=amplifyProp, min_length=minCnvLength, max_length=maxCnvLength, exp_length=expCnvLength, append=False)
        log("Simulated CNV list is at " + repr(cnvList))

    log("Start simulating diploid genomes")
    GenerateGenomes(genome_file=genomeFile, cnv_list_file=cnvList, paternal_file=paternalGenomeFile, maternal_file=maternalGenomeFile)
    
    log("Finished simulating genomes, starting simulating NGS reads")
    CallART(genome_file=paternalGenomeFile, output_file=paternalReadsFile, read_length=readLength, fold_coverage=coverage)
    CallART(genome_file=maternalGenomeFile, output_file=maternalReadsFile, read_length=readLength, fold_coverage=coverage)
    
    log("Finished simulating NGS reads, starting merging reads from paternal and maternal genomes")
    paternalFile1 = paternalReadsFile + '1.fq'
    paternalFile2 = paternalReadsFile + '2.fq'
    maternalFile1 = maternalReadsFile + '1.fq'
    maternalFile2 = maternalReadsFile + '2.fq'
    outReadsFile1 = os.path.join(outDir, prefix+'_1.fq')
    outReadsFile2 = os.path.join(outDir, prefix+'_2.fq')
    MergeReads(paternalFile1, maternalFile1, outReadsFile1)
    MergeReads(paternalFile2, maternalFile2, outReadsFile2)

    log("Start compressing merged fastq files")
    subprocess.call(['gzip', '--force', outReadsFile1])
    subprocess.call(['gzip', '--force', outReadsFile2])

    # remove unwanted intermediate files
    subprocess.call(['rm', paternalFile1, paternalFile2, maternalFile1, maternalFile2])
    log("All task finished!")
    

if __name__ == '__main__':
    main()