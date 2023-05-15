#! /usr/bin/env python

import sys
from collections import defaultdict
import random
import numpy as np
import copy
import subprocess
import datetime
import csv

def log(message):
    print('[CNV Simulator {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + "] " + message)


def overlapCNV(cnv_region, cnv_list):
    """Check whether the single CNV overlaps with a list of CNVs"""
    if not cnv_list:
        return False
    c1, s1, e1 = cnv_region
    for x in cnv_list:
        c2, s2, e2 = x[0], int(x[1]), int(x[2])
        if c1 == c2 and (s1 < s2 < e1 or s1 < e2 < e1):
            return True
        
    return False


def SimulateCNV(chrom, access_region_by_chr, cnv_num, p_amplify=0.5, min_length=5000, max_length=2000000, exp_length=200000, append=False, cnv_list=None):
    """
    Simulate CNVs in accessible regions of single chromosomes
    - access_region_by_chr: a list of accessible regions for a single chromosome
    - cnv_num: the number of CNVs need to be simulated in that chromosome
    return: a list with elements as: [start, end, CN]
    """

    cnvRegions = []

    # if the file of CNV list is provided, make sure newly generated CNVs do not overlap with them
    if append:
        tmpCnvRegions = cnv_list
    else:
        tmpCnvRegions = []

    while cnv_num:
        # randomly choose a region to generate CNV
        accessRegion = access_region_by_chr[random.randrange(len(access_region_by_chr))]
        regionStart, regionEnd = accessRegion
        # randomly generate the length of the CNV, follow exponential distribution, scale=200,000
        tmpCnvLength = round(np.random.exponential(exp_length))
        cnvLength = tmpCnvLength + min_length
        cnvLength = max_length if cnvLength > max_length else cnvLength
        # randomly generate the start location of the CNV
        startLoc = random.randint(regionStart, regionEnd)
        
        # if the CNV exceeds the boundary of the accessible region, simulate again
        if startLoc + cnvLength > regionEnd:
            continue
        cnvRegion = [chrom, startLoc, startLoc + cnvLength]

        # make sure newly produced CNVs not overlap with old CNVs
        if not overlapCNV(cnvRegion, tmpCnvRegions):
            tmpCnvRegions.append(cnvRegion)
        else:
            continue

        # simulate the copy number for the CNV
        if random.random() < p_amplify:
            cn = random.choices([0, 1], [1, 2], k=1)[0]
        else:
            cn = random.choices([3, 4, 5], [4, 4, 2], k=1)[0]
        
        cnvRegion.append(cn)
        cnvRegions.append(cnvRegion)
        cnv_num -= 1

    return sorted(cnvRegions)


def GenerateCNVList(access_file, output_file, cnv_num, p_amplify=0.5, min_length=5000, max_length=500000, exp_length=200000, append=False):
    """
    Randomly generate a CNV list for next step to produce new genomes with CNVs
    - access_file: includes genomic accessible regions
    - output_file: output file with a list of CNVs
    - cnv_num: total number of CNVs
    - p_amplify: proportion of copy number amplifications
    - min_length: minimum length of CNV
    - max_length: maximum length of CNV
    - exp_length: lambda of exponential distribution to simulate the length of CNV
    - append: whether append newly simualted CNVs into existing file
    return: None, write a CNV list to a file with columns as: chr, start, end, CN
    """

    # get the total length of chromosomal accessible regions
    chrs = ['chr'+str(i) for i in range(1,23)] + ['chrX']
    chrLenDict = {chrom:0 for chrom in chrs}
    accessRegionByChr = defaultdict(list)
    with open(access_file) as f:
        for line in f:
            accessRegion = line.strip().split('\t')
            chrom = accessRegion[0]
            if chrom not in chrs:
                continue
            start, end = int(accessRegion[1]), int(accessRegion[2])
            accessRegionByChr[chrom].append([start, end])
            chrLenDict[chrom] += end - start
    
    # estimate how many CNVs shoul be simulated in each chromosome according to total CNV number 
    # and chromosome length
    cnvNumDict = {}
    for i in range(1,len(chrs)):
        chrom = chrs[i]
        cnvNumDict[chrom] = round(chrLenDict[chrom] * cnv_num / sum(chrLenDict.values()))
    
    cnvNumDict['chr1'] = cnv_num - sum(cnvNumDict.values())

    # simulate CNVs by chromosome, depending on the mode
    if append:
        mode = 'a'
        cnv_list = list(csv.reader(open(output_file, 'r'), delimiter='\t'))
    else:
        mode = 'w'
        cnv_list = None
        
    simuCnvByChr = {}
    for chrom in cnvNumDict.keys():
        cnvNum = cnvNumDict[chrom]
        accessRegion = accessRegionByChr[chrom]
        tmpCnvList = SimulateCNV(chrom, accessRegion, cnvNum, p_amplify, min_length, max_length, exp_length, append=append, cnv_list=cnv_list)
        simuCnvByChr[chrom] = tmpCnvList
    
    # get final simulated CNV list
    finalCnvList = [cnvs for chrom in chrs for cnvs in simuCnvByChr[chrom]]
    
    # write or append CNV list into a file
    csv.writer(open(output_file, mode=mode, newline=''), delimiter='\t').writerows(finalCnvList)
        

def GenerateGenomes(genome_file, cnv_list_file, paternal_file='paternal.fa', maternal_file='maternal.fa'):
    """
    - genome_file: a human genome fasta file
    - cnv_list_file: a file with a list of CNV regions
    return: None, write two genomes (paternal and maternal) with varied sequences into files
    """

    chrs = ['chr'+str(i) for i in range(1,23)] + ['chrX']
    
    # read CNV list file
    cnv_list = list(csv.reader(open(cnv_list_file, 'r'), delimiter='\t'))

    # determine the CNV allele by chromosome, to make it simple, put all variations into one chromosome    
    # also aggregate CNV list by chromosome
    cnvListByChr = defaultdict(list)
    paternalAllele, maternalAllele = defaultdict(list), defaultdict(list)
    for x in cnv_list:
        
        tmpChrom, cn = x[0], int(x[3])
        cnvListByChr[tmpChrom].append(x[1:3])
        if cn == 0:
            paternalAllele[tmpChrom].append(0)
            maternalAllele[tmpChrom].append(0)
        else:
            paternalAllele[tmpChrom].append(1)
            maternalAllele[tmpChrom].append(cn-1)

    # reverse the CNV in the list for easier sequence manipulation
    for x in chrs:
        paternalAllele[x] = paternalAllele[x][::-1]
        maternalAllele[x] = maternalAllele[x][::-1]
    
    # print('paternal CNV allele is: ', paternalAllele)
    # print('maternal CNV allele is: ', maternalAllele)

    # read chromosomal sequences into a dictionary. very time-consuming
    log("Start loading genome")
    chrSeq = {}
    sequence = ''
    chrom = 'none'
    with open(genome_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if line[1:] not in chrs:
                    continue
                else:
                    chrSeq[chrom] = sequence
                    chrom = line[1:]
                    sequence = ''
            else:
                sequence += line
        
        chrSeq[chrom] = sequence
    
    del chrSeq['none']

    # modify the chromosomal sequences by CNVs
    out_paternal_file = open(paternal_file, 'w')
    out_maternal_file = open(maternal_file, 'w')
    for chrom in chrs:
        log("Adding CNVs to chromosome " + chrom)
        tmpSeq = chrSeq[chrom]
        tmpSeqListPaternal = list(tmpSeq)
        tmpSeqListMaternal = copy.deepcopy(tmpSeqListPaternal)
        for j, cnvRegion in enumerate(cnvListByChr[chrom][::-1]):
            start, end = int(cnvRegion[0]), int(cnvRegion[1])
            # print('Start is {:d}, end is {:d}'.format(start, end))
            paternalcn = paternalAllele[chrom][j]
            maternalcn = maternalAllele[chrom][j]
            if paternalcn != 1:
                del tmpSeqListPaternal[start:end]
            if maternalcn > 1:
                tmpSeqListMaternal = tmpSeqListMaternal[:start] + tmpSeqListMaternal[start:end]*maternalcn + tmpSeqListMaternal[end:]
            elif maternalcn < 1:
                del tmpSeqListMaternal[start:end]
            else:
                pass
        
        print('>'+chrom, file=out_paternal_file)
        print(''.join(tmpSeqListPaternal), file=out_paternal_file)
        print('>'+chrom, file=out_maternal_file)
        print(''.join(tmpSeqListMaternal), file=out_maternal_file)


def CallART(genome_file, output_file, read_length, fold_coverage=1):
    """
    Call ART to generate artificial reads
    - genome_file: reference genome file in FASTA format
    - output_file: output file name
    - read_length: the read length
    - fold_coverage: fold coverage for the reads
    return: None
    """
    # make sure art_illumina is in global environment
    subprocess.call(["art_illumina", \
                     "-na", \
                     "-ss", "HSXn", \
                     "-i", genome_file, \
                     "-p", \
                     "-m", "800", \
                     "-s", "10", \
                     "-l", str(read_length), \
                     "-f", str(fold_coverage), \
                     "-o", output_file], stderr=None)
    

def MergeReads(infile1, infile2, outfile):
    """Concatenate two files by Linux Shell cat"""
    with open(outfile, 'w') as f:
        subprocess.call(['cat', infile1, infile2], stdout=f)


if __name__ == "__main__":

    accessFile = sys.argv[1]
    genomeFile = sys.argv[2]
    outputFile = sys.argv[3]
    paternalFile = sys.argv[4]
    maternalFile = sys.argv[5]

    GenerateCNVList(access_file=accessFile, output_file=outputFile, cnv_num=100)
    GenerateGenomes(genomeFile, outputFile, paternal_file=paternalFile, maternal_file=maternalFile)
