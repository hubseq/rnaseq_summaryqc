import os, subprocess, sys
import pandas as pd
import numpy as np
import zipfile  # requires Python > 3.2
# sys.path.append('global_utils/src/')
sys.path.append('/global_utils/src/')
import module_utils
import aws_s3_utils
import file_utils


def parse_fastqc_fastqc(fastqc_file):
    """ Specifically parse FASTQC data file

        {'per_base_quality': {'mean': [...], 'median': [...], 'pos': [...], 'lowerquartile': [...], 'upperquartile': [...]},
         'per_tile_quality': {'tile': [...], 'base': [...], ''mean': [...]},
        }
    """
    fastqc_stats = {'per_base_quality': {'mean': [], 'median': [], 'pos': [], 'lowerquartile': [], 'upperquartile': []},
                    'per_tile_quality': {'tile': [], 'base': [], 'mean': []},
                    'per_base_content': {'pos': [], 'A': [], 'T': [], 'G': [], 'C': []},
                    'per_base_GC_content': {'pos': [], 'count': [], 'total_count': 0},
                    'per_base_N_content': {'pos': [], 'perc': []},
                    'deduplicated_perc': 100.0,
                    'adapter_content': {'pos': [], 'illumina_universal_adapter': [], 'illumina_small_rna_3p_adapter': [],
                                        'illumina_small_rna_5p_adapter': [], 'illumina_nextera_transposase': [],
                                        'solid_small_rna_adapter': []},
                    'read_length': 0,
                    'total_raw_reads': 0
                    }
    with open(fastqc_file,'r') as f:
        while True:
            r = f.readline()
            if r=='':
                break
            else:
                if r.startswith('>>Per base sequence quality'):
                    r = f.readline()  #Base  Mean  Median  Lower Quartile  Upper Quartile  10th Percentile 90th Percentile
                    r = f.readline()  # first position
                    while not r.startswith('>>END_MODULE'):
                        if not r.isspace():
                            rt = r.rstrip(' \t\n').split('\t')
                            fastqc_stats['per_base_quality']['pos'].append(int(rt[0]))
                            fastqc_stats['per_base_quality']['mean'].append(float(rt[1]))
                            fastqc_stats['per_base_quality']['median'].append(float(rt[2]))
                            fastqc_stats['per_base_quality']['lowerquartile'].append(float(rt[3]))
                            fastqc_stats['per_base_quality']['upperquartile'].append(float(rt[4]))
                        r = f.readline()
                elif r.startswith('Sequence length'):
                    rt = r.rstrip(' \t\n').split('\t')
                    fastqc_stats['read_length'] = int(rt[1])
                elif r.startswith('Total Sequences'):
                    rt = r.rstrip(' \t\n').split('\t')
                    fastqc_stats['total_raw_reads'] = int(rt[1])                    
                elif r.startswith('>>Per tile sequence quality'):
                    r = f.readline()  #Tile  Base  Mean
                    r = f.readline()  # first position
                    while not r.startswith('>>END_MODULE'):
                        if not r.isspace():
                            rt = r.rstrip(' \t\n').split('\t')
                            fastqc_stats['per_tile_quality']['tile'].append(int(rt[0]))
                            fastqc_stats['per_tile_quality']['base'].append(int(rt[1]))
                            fastqc_stats['per_tile_quality']['mean'].append(float(rt[2]))
                        r = f.readline()                            
                elif r.startswith('>>Per base sequence content'):
                    r = f.readline()  #Base   G       A       T       C
                    r = f.readline()  # first position
                    while not r.startswith('>>END_MODULE'):
                        if not r.isspace():
                            rt = r.rstrip(' \t\n').split('\t')
                            fastqc_stats['per_base_content']['pos'].append(int(rt[0]))
                            fastqc_stats['per_base_content']['G'].append(float(rt[1]))
                            fastqc_stats['per_base_content']['A'].append(float(rt[2]))
                            fastqc_stats['per_base_content']['T'].append(float(rt[3]))
                            fastqc_stats['per_base_content']['C'].append(float(rt[4]))
                        r = f.readline()
                elif r.startswith('>>Per sequence GC content'):
                    r = f.readline()  #GC Content     Count
                    r = f.readline()  # first position
                    total_count = 0
                    while not r.startswith('>>END_MODULE'):
                        if not r.isspace():
                            rt = r.rstrip(' \t\n').split('\t')
                            fastqc_stats['per_base_GC_content']['pos'].append(int(rt[0]))
                            fastqc_stats['per_base_GC_content']['count'].append(float(rt[1]))
                            total_count += int(float(rt[1]))
                        r = f.readline()                        
                    fastqc_stats['per_base_GC_content']['total_count'] = total_count
                elif r.startswith('>>Per base N content'):
                    r = f.readline()  #Base   N-Count (%)
                    r = f.readline()  # first position
                    while not r.startswith('>>END_MODULE'):
                        if not r.isspace():
                            rt = r.rstrip(' \t\n').split('\t')
                            fastqc_stats['per_base_N_content']['pos'].append(int(rt[0]))
                            fastqc_stats['per_base_N_content']['perc'].append(float(rt[1]))
                        r = f.readline()
                elif r.startswith('>>Sequence Duplication Levels'):
                    r = f.readline()  # total deduplicated %
                    rt = r.rstrip(' \t\n').split('\t')
                    fastqc_stats['deduplicated_perc'] = float(rt[1])
                elif r.startswith('>>Adapter content'):
                    r = f.readline()  #Position       Illumina Universal Adapter      Illumina Small RNA 3' Adapter   Illumina Small RNA 5' Adapter   Nextera Transposase Sequence    SOLID Small RNA Adapter
                    r = f.readline()  # first position
                    while not r.startswith('>>END_MODULE'):
                        if not r.isspace():
                            rt = r.rstrip(' \t\n').split('\t')
                            fastqc_stats['adapter_content']['pos'].append(int(rt[0]))
                            fastqc_stats['adapter_content']['illumina_universal_adapter'].append(rt[0])
                            fastqc_stats['adapter_content']['illumina_small_rna_3p_adapter'].append(rt[1])
                            fastqc_stats['adapter_content']['illumina_small_rna_5p_adapter'].append(rt[2])
                            fastqc_stats['adapter_content']['illumina_nextera_transposase'].append(rt[3])
                            fastqc_stats['adapter_content']['solid_small_rna_adapter'].append(rt[4])
                        r = f.readline()
    return fastqc_stats


def summary_fastqc_fastqc(fastqc_stats):
    """ After getting FASTQC raw stats, create a JSON that summarizes information.

    fastqc_stats = {'per_base_quality': {'mean': [], 'median': [], 'pos': [], 'lowerquartile': [], 'upperquartile': []},
                    'per_tile_quality': {'tile': [], 'base': [], 'mean': []},
                    'per_base_content': {'A': [], 'T': [], 'G': [], 'C': []},
                    'per_base_GC_content': {'pos': [], 'count': [], 'total_count': 0},
                    'per_base_N_content': {'pos': [], 'perc': []},
                    'deduplicated_perc': 100.0,
                    'adapter_content': {'pos': [], 'illumina_universal_adapter': [], 'illumina_small_rna_3p_adapter': [],
                                        'illumina_small_rna_5p_adapter': [], 'illumina_nextera_transposase': [],
                                        'solid_small_rna_adapter': []}
                    }
    """
    fastqc_summary = {'per_base_quality': {'result': 'PASS', 'failed_pos_low': [], 'suggest': [], 'value': [], 'threshold': 'Q>25' },
                      'per_base_content': {'result': 'PASS', 'failed_pos': [], 'failed_pos_low': [], 'failed_pos_high': [], 'failed_pos_low_nt': [], 'failed_pos_high_nt': [], 'suggest': [], 'value': [], 'threshold': '15-35% per base' },
                      'per_base_GC_content': {'result': 'PASS', 'suggest': [], 'value': [], 'threshold': '30-70%' },
                      'per_base_N_content': {'result': 'PASS', 'failed_pos': [], 'suggest': [], 'value': [], 'threshold': '<10%' },
                      'deduplicated_perc': {'result': 'PASS', 'suggest': [], 'value': [], 'threshold': '>20%' },
                      'adapter_content': {'result': 'PASS', 'failed_pos': [], 'suggest': [], 'value': [], 'threshold': '<10%' }
                      }
    
    ## report any mean base quality < 25
    for i in range(0,len(fastqc_stats['per_base_quality']['pos'])):
        if float(fastqc_stats['per_base_quality']['mean'][i]) < 25:
            fastqc_summary['per_base_quality']['failed_pos_low'].append(int(fastqc_stats['per_base_quality']['pos'][i]))
    if len(fastqc_summary['per_base_quality']['failed_pos_low']) > 0:
        fastqc_summary['per_base_quality']['result'] = 'WARNING'
        
        current_pos = 0
        for i in range(0,len(fastqc_summary['per_base_quality']['failed_pos_low'])):
            if fastqc_summary['per_base_quality']['failed_pos_low'][i] == current_pos + 1:
                current_pos += 1
            else:
                break
        if current_pos > 0:
            fastqc_summary['per_base_quality']['suggest'].append("Read Sequencing - low sequencing quality detected on the 5'end. Suggest 5'end read trimming.")
            fastqc_summary['per_base_quality']['value'].append(current_pos)

        pos_len = len(fastqc_summary['per_base_quality']['failed_pos_low'])
        max_pos = max(fastqc_stats['per_base_quality']['pos'])
        current_pos = max_pos
        for i in range(0,pos_len):
            if fastqc_summary['per_base_quality']['failed_pos_low'][pos_len-i-1] == current_pos:
                current_pos = current_pos - 1
            else:
                break
        if current_pos < max_pos:
            fastqc_summary['per_base_quality']['suggest'].append("Read Sequencing - low sequencing quality detected on the 3'end. Suggest 3'end read trimming.")
            fastqc_summary['per_base_quality']['value'].append(max_pos-current_pos)
            
    ## report any per-base content nt bias > 35% or < 15%
    for i in range(0,len(fastqc_stats['per_base_content']['pos'])):
        if float(fastqc_stats['per_base_content']['A'][i]) > 35.0:
            fastqc_summary['per_base_content']['failed_pos_high'].append(fastqc_stats['per_base_content']['pos'][i])
            fastqc_summary['per_base_content']['failed_pos_high_nt'].append('A')
            fastqc_summary['per_base_content']['failed_pos'].append(fastqc_stats['per_base_content']['pos'][i])
        elif float(fastqc_stats['per_base_content']['A'][i]) < 15.0:
            fastqc_summary['per_base_content']['failed_pos_low'].append(fastqc_stats['per_base_content']['pos'][i])
            fastqc_summary['per_base_content']['failed_pos_low_nt'].append('A')
            fastqc_summary['per_base_content']['failed_pos'].append(fastqc_stats['per_base_content']['pos'][i])
            
        if float(fastqc_stats['per_base_content']['T'][i]) > 35.0:
            fastqc_summary['per_base_content']['failed_pos_high'].append(fastqc_stats['per_base_content']['pos'][i])
            fastqc_summary['per_base_content']['failed_pos_high_nt'].append('T')
            fastqc_summary['per_base_content']['failed_pos'].append(fastqc_stats['per_base_content']['pos'][i])
        elif float(fastqc_stats['per_base_content']['T'][i]) < 15.0:
            fastqc_summary['per_base_content']['failed_pos_low'].append(fastqc_stats['per_base_content']['pos'][i])
            fastqc_summary['per_base_content']['failed_pos_low_nt'].append('T')
            fastqc_summary['per_base_content']['failed_pos'].append(fastqc_stats['per_base_content']['pos'][i])
        
        if float(fastqc_stats['per_base_content']['C'][i]) > 35.0:
            fastqc_summary['per_base_content']['failed_pos_high'].append(fastqc_stats['per_base_content']['pos'][i])
            fastqc_summary['per_base_content']['failed_pos_high_nt'].append('C')
            fastqc_summary['per_base_content']['failed_pos'].append(fastqc_stats['per_base_content']['pos'][i])
        elif float(fastqc_stats['per_base_content']['C'][i]) < 15.0:
            fastqc_summary['per_base_content']['failed_pos_low'].append(fastqc_stats['per_base_content']['pos'][i])
            fastqc_summary['per_base_content']['failed_pos_low_nt'].append('C')
            fastqc_summary['per_base_content']['failed_pos'].append(fastqc_stats['per_base_content']['pos'][i])

        if float(fastqc_stats['per_base_content']['G'][i]) > 35.0:
            fastqc_summary['per_base_content']['failed_pos_high'].append(fastqc_stats['per_base_content']['pos'][i])
            fastqc_summary['per_base_content']['failed_pos_high_nt'].append('G')
            fastqc_summary['per_base_content']['failed_pos'].append(fastqc_stats['per_base_content']['pos'][i])
        elif float(fastqc_stats['per_base_content']['G'][i]) < 15.0:
            fastqc_summary['per_base_content']['failed_pos_low'].append(fastqc_stats['per_base_content']['pos'][i])
            fastqc_summary['per_base_content']['failed_pos_low_nt'].append('G')
            fastqc_summary['per_base_content']['failed_pos'].append(fastqc_stats['per_base_content']['pos'][i])
    
    # remove repeat positions, preserving order
    fastqc_summary['per_base_content']['failed_pos'] = list(set(fastqc_summary['per_base_content']['failed_pos']))

    if len(fastqc_summary['per_base_content']['failed_pos']) > 0:
        fastqc_summary['per_base_content']['result'] = 'WARNING'
        
        current_pos = 0
        for i in range(0,len(fastqc_summary['per_base_content']['failed_pos'])):
            if fastqc_summary['per_base_content']['failed_pos'][i] == current_pos + 1:
                current_pos += 1
            else:
                break
        if current_pos > 0:
            fastqc_summary['per_base_content']['suggest'].append("Read sequencing - detected nucleotide bias on the 5'end. Suggest 5'end read trimming.")
            fastqc_summary['per_base_content']['value'].append(current_pos)                             
        
        pos_len = len(fastqc_summary['per_base_content']['failed_pos_low'])
        max_pos = max(fastqc_stats['per_base_content']['pos'])
        current_pos = max_pos
        for i in range(0,pos_len):
            if fastqc_summary['per_base_content']['failed_pos_low'][pos_len-i-1] == current_pos:
                current_pos -=1
            else:
                break
        if current_pos < max_pos:
            fastqc_summary['per_base_content']['suggest'].append("Read sequencing - detected nucleotide bias on the 3'end. Suggest 3'end read trimming.")
            fastqc_summary['per_base_content']['value'].append(max_pos-current_pos)

    ## report extreme GC% bias
    max_GC_perc_index = fastqc_stats['per_base_GC_content']['count'].index(max(fastqc_stats['per_base_GC_content']['count']))
    max_GC_perc = fastqc_stats['per_base_GC_content']['pos'][max_GC_perc_index]
    if max_GC_perc < 20 or max_GC_perc > 80:
        fastqc_summary['per_base_GC_content']['result'] = 'FAIL'
        fastqc_summary['per_base_GC_content']['suggest'] = 'Read sequencing - detected extreme GC bias.'
        fastqc_summary['per_base_GC_content']['value'] = max_GC_perc
    elif max_GC_perc < 30 or max_GC_perc > 70:
        fastqc_summary['per_base_GC_content']['result'] = 'WARNING'
        fastqc_summary['per_base_GC_content']['suggest'] = 'Read sequencing - detected moderate GC bias.'
        fastqc_summary['per_base_GC_content']['value'] = max_GC_perc

    ## report any positions that have a high number of Ns
    max_pos = max(fastqc_stats['per_base_N_content']['pos'])    
    for i in range(0,len(fastqc_stats['per_base_N_content']['pos'])):
        if float(fastqc_stats['per_base_N_content']['perc'][i]) > 10:
            fastqc_summary['per_base_N_content']['failed_pos'].append(int(fastqc_stats['per_base_N_content']['pos'][i]))
    suggested_mismatch_rate = int(100.0*len(fastqc_summary['per_base_N_content']['failed_pos'])/max_pos)
    if len(fastqc_summary['per_base_N_content']['failed_pos']) > 0 and suggested_mismatch_rate <= 15:
        fastqc_summary['per_base_N_content']['result'] = 'WARNING'
        fastqc_summary['per_base_N_content']['suggest'] = 'Read sequencing - detected greater than 10% unknown base calls (N) at certain base positions. Suggest increasing allowed mismatch rate during alignment.'
        fastqc_summary['per_base_N_content']['value'].append(fastqc_summary['per_base_N_content']['failed_pos'])
        fastqc_summary['per_base_N_content']['value'].append(suggested_mismatch_rate)
    elif len(fastqc_summary['per_base_N_content']['failed_pos']) > 0 and suggested_mismatch_rate > 15:
        fastqc_summary['per_base_N_content']['result'] = 'FAIL'
        fastqc_summary['per_base_N_content']['suggest'] = 'Read sequencing - detected greater than 10% unknown base calls (N) at too many positions. Alignment would require an allowed mismatch rate of greater than 15%.'
        fastqc_summary['per_base_N_content']['value'].append(fastqc_summary['per_base_N_content']['failed_pos'])
        fastqc_summary['per_base_N_content']['value'].append(suggested_mismatch_rate)
    
    ## report any deduplication less than 20%
    if float(fastqc_stats['deduplicated_perc']) < 20.0:
        fastqc_summary['deduplicated_perc']['result'] = 'WARNING'
        fastqc_summary['deduplicated_perc']['suggest'] = 'Read sequencing - only {}% of reads will remain after removing read duplicates (reads that will map to same genomic coordinates). This suggests low sequencing complexity or low concentration of sequenced sample.'
    
    ## report any required adapter trimming - 'adapter_content'
    for i in range(0,len(fastqc_stats['adapter_content']['pos'])):
        if float(fastqc_stats['adapter_content']['illumina_universal_adapter'][i]) > 20.0 or \
           float(fastqc_stats['adapter_content']['illumina_small_rna_3p_adapter'][i]) > 20.0 or \
           float(fastqc_stats['adapter_content']['illumina_small_rna_5p_adapter'][i]) > 20.0 or \
           float(fastqc_stats['adapter_content']['illumina_nextera_transposase'][i]) > 20.0 or \
           float(fastqc_stats['adapter_content']['solid_small_rna_adapter'][i]) > 20.0:
            fastqc_summary['adapter_content']['failed_pos'].append(fastqc_stats['adapter_content']['pos'])
    
    if len(fastqc_summary['adapter_content']['failed_pos']) > 0:
        fastqc_summary['adapter_content']['result'] = 'WARNING'
        
        current_pos = 0
        for i in range(0,len(fastqc_summary['adapter_content']['failed_pos'])):
            if fastqc_summary['adapter_content']['failed_pos'][i] == current_pos + 1:
                current_pos += 1
            else:
                break
        if current_pos > 0:
            fastqc_summary['adapter_content']['suggest'].append("Read sequencing - detected sequencing adapter on the 5'end. Suggest 5'end read trimming.")
            fastqc_summary['adapter_content']['value'] = current_pos
                                                            

        pos_len = len(fastqc_summary['adapter_content']['failed_pos'])
        max_pos = max(fastqc_stats['adapter_content']['pos'])
        current_pos = max_pos
        for i in range(0,pos_len):
            if fastqc_summary['adapter_content']['failed_pos_low'][pos_len-i-1] == current_pos:
                current_pos -=1
            else:
                break
        if current_pos < max_pos:
            fastqc_summary['adapter_content']['suggest'].append("Read sequencing - detected sequencing adapter on the 3'end. Suggest 3'end read trimming.")
            fastqc_summary['adapter_content']['value'].append(max_pos-current_pos)
    
    return fastqc_summary


def parse_aligner_rnastar( aligner_file ):
    """ Specifically parse RNA-STAR aligner data file

    """
    aligner_stats = {'input_reads': 0, 'percent_mapped': 0, 'percent_uniquely_mapped': 0, 'percent_mismatch_rate': 0, 
                     'percent_deletion_rate': 0, 'percent_insertion_rate': 0, 'average_deletion_length': 0, 'average_insertion_length': 0
                    }

    percent_mapped = 0
    
    with open(aligner_file,'r') as f:
        while True:
            r = f.readline()
            if r=='':
                break
            else:
                if 'Number of input reads' in r:
                    aligner_stats['input_reads'] = int(r.rstrip(' \t\n').split('\t')[-1])
                elif 'Uniquely mapped reads %' in r:
                    unique_mapped = float(r.rstrip(' \t\n').split('\t')[-1].rstrip('%'))
                    aligner_stats['percent_uniquely_mapped'] = unique_mapped
                    percent_mapped += unique_mapped
                elif '% of reads mapped to multiple loci' in r:
                    multiple_mapped = float(r.rstrip(' \t\n').split('\t')[-1].rstrip('%'))
                    percent_mapped += multiple_mapped
                elif '% of reads mapped to too many loci' in r:
                    too_many_mapped = float(r.rstrip(' \t\n').split('\t')[-1].rstrip('%'))
                    percent_mapped += too_many_mapped
                elif 'Mismatch rate per base, %' in r:
                    aligner_stats['percent_mismatch_rate'] = float(r.rstrip(' \t\n').split('\t')[-1].rstrip('%'))
                elif 'Deletion rate per base' in r:
                    aligner_stats['percent_deletion_rate'] = float(r.rstrip(' \t\n').split('\t')[-1].rstrip('%'))
                elif 'Insertion rate per base' in r:
                    aligner_stats['percent_insertion_rate'] = float(r.rstrip(' \t\n').split('\t')[-1].rstrip('%'))
                elif 'Deletion average length' in r:
                    aligner_stats['average_deletion_length'] = float(r.rstrip(' \t\n').split('\t')[-1])
                elif 'Insertion average length' in r:
                    aligner_stats['average_insertion_length'] = float(r.rstrip(' \t\n').split('\t')[-1])
        aligner_stats['percent_mapped'] = percent_mapped
    return aligner_stats

def summary_aligner_rnastar( aligner_stats ):
    """ After getting aligner stats, generate summary

        aligner_stats = {'input_reads': 0, 'percent_mapped': 0, 'percent_uniquely_mapped': 0, 'percent_mismatch_rate': 0, 
                     'percent_deletion_rate': 0, 'percent_insertion_rate': 0, 'average_deletion_length': 0, 'average_insertion_length': 0
                    }
    """
    aligner_summary = {'percent_mapped': {'result': 'PASS', 'suggest': '', 'value': 'NA', 'threshold': '>50%' },
                       'percent_uniquely_mapped': {'result': 'PASS', 'suggest': '', 'value': 'NA', 'threshold': '>50% of mapped' },
                       'percent_mismatch_rate': {'result': 'PASS', 'suggest': '', 'value': 'NA', 'threshold': '<15%' },
                       'percent_deletion_rate': {'result': 'PASS', 'suggest': '', 'value': 'NA', 'threshold': '<10%' },
                       'percent_insertion_rate': {'result': 'PASS', 'suggest': '', 'value': 'NA', 'threshold': '<10%' }
                      }

    print('ALIGNER STATS: '+str(aligner_stats))
    
    ## report if % mapped < 50%
    if float(aligner_stats['percent_mapped']) < 50.0:
        aligner_summary['percent_mapped']['result'] = 'WARNING'
        aligner_summary['percent_mapped']['suggest'] = 'Alignment mapping rate is less than 50%. Consider read trimming or increasing allowed mismatch rate as suggested through FASTQ QC, or ignoring this sample if possible.'
        aligner_summary['percent_mapped']['value'] = round(float(aligner_stats['percent_mapped']),2)
    elif float(aligner_stats['percent_mapped']) < 20.0:
        aligner_summary['percent_mapped']['result'] = 'FAIL'
        aligner_summary['percent_mapped']['suggest'] = 'Alignment mapping rate is less than 20%. Consider removing this sample from further analysis.'
        aligner_summary['percent_mapped']['value'] = round(float(aligner_stats['percent_mapped']),2)        

    ## report if % uniquely mapped is < 50% of mapped
    if float(aligner_stats['percent_mapped']) > 0 and float(aligner_stats['percent_uniquely_mapped'])/float(aligner_stats['percent_mapped']) < 0.5:
        aligner_summary['percent_uniquely_mapped']['result'] = 'WARNING'
        aligner_summary['percent_uniquely_mapped']['suggest'] = 'More than 50% of mapped reads are aligning to multiple regions. If this is an issue, consider more stringent mapping parameters.'
        aligner_summary['percent_uniquely_mapped']['value'] = float(aligner_stats['percent_uniquely_mapped'])/float(aligner_stats['percent_mapped']) if float(aligner_stats['percent_mapped']) > 0 else 0

    ## report if mismatch rate > 15%
    if float(aligner_stats['percent_mismatch_rate']) > 15.0:
        aligner_summary['percent_mismatch_rate']['result'] = 'WARNING'
        aligner_summary['percent_mismatch_rate']['suggest'] = 'Mismatch rate > 15%. Consider more stringent mapping parameters, or re-sequencing if possible.'
        aligner_summary['percent_mismatch_rate']['value'] = round(float(aligner_stats['percent_mismatch_rate']),2)

    ## report if deletion rate > 10%
    if float(aligner_stats['percent_deletion_rate']) > 10.0:
        aligner_summary['percent_deletion_rate']['result'] = 'WARNING'
        aligner_summary['percent_deletion_rate']['suggest'] = 'Deletion rate > 10%. Consider more stringent mapping parameters, or re-sequencing if possible.'
        aligner_summary['percent_deletion_rate']['value'] = round(float(aligner_stats['percent_deletion_rate']),2)

    ## report if insertion rate > 10%
    if float(aligner_stats['percent_insertion_rate']) > 10.0:
        aligner_summary['percent_insertion_rate']['result'] = 'WARNING'
        aligner_summary['percent_insertion_rate']['suggest'] = 'Insertion rate > 10%. Consider more stringent mapping parameters, or re-sequencing if possible.'
        aligner_summary['percent_insertion_rate']['value'] = round(float(aligner_stats['percent_insertion_rate']),2)
    
    return aligner_summary


def parse_alignqc_rnaseqc( alignqc_file ):
    """ Specifically parse RNA-STAR aligner data file
    
    """
    alignqc_stats = {'sample_file': '', 'mapping_rate_in_file': 0, 'unique_mapping_rate_of_mapped_in_file': 0, 'duplicate_mapping_rate_of_mapped_in_file': 0,
                     'base_mismatch_rate': 0, 'R1_mapping_rate': 0, 'R2_mapping_rate': 0, 'R1_mismatch_rate': 0, 'R2_mismatch_rate': 0,
                     'exon_mapping_rate': 0, 'intron_mapping_rate': 0, 'intergenic_mapping_rate': 0, 'intragenic_mapping_rate': 0,
                     'ambiguous_alignment_rate': 0, 'rRNA_mapping_rate': 0, 'high_qualiity_rate': 0,
                     'R1_sense_strand_mapping_rate': 0, 'R2_sense_strand_mapping_rate': 0,
                     'average_insertions_deletions_gaps_per_read': 0, 'failed_vendor_QC': 0, 'total_reads': 0, 'percent_failed_vendor_QC': 0,
                     'low_mapping_quality': 0, 'percent_low_mapping_quality': 0, 'genes_detected': 0,
                     'mean_gene_3p_bias_in_read_location': 0, 'median_gene_3p_bias_in_read_location': 0,
                     'median_transcript_gene_coverage': 0 }
    
    with open(alignqc_file,'r') as f:
        while True:
            r = f.readline()
            if r=='':
                break
            elif not r.isspace():
                rt = r.rstrip(' \t\n').split('\t')
                if len(rt) > 1:
                    if rt[0]=='Sample':
                        alignqc_stats['sample_file'] = rt[1]
                    elif rt[0]=='Mapping Rate':
                        alignqc_stats['mapping_rate_in_file'] = float(rt[1])
                    elif rt[0]=='Unique Rate of Mapped':
                        alignqc_stats['unique_mapping_rate_of_mapped_in_file'] = float(rt[1])
                    elif rt[0]=='Duplicate Rate of Mapped':
                        alignqc_stats['duplicate_mapping_rate_of_mapped_in_file'] = float(rt[1])
                    elif rt[0]=='Base Mismatch':
                        alignqc_stats['base_mismatch_rate'] = float(rt[1])
                    elif rt[0]=='End 1 Mapping Rate':
                        alignqc_stats['R1_mapping_rate'] = float(rt[1])
                    elif rt[0]=='End 2 Mapping Rate':
                        alignqc_stats['R2_mapping_rate'] = float(rt[1])
                    elif rt[0]=='End 1 Mismatch Rate':
                        alignqc_stats['R1_mismatch_rate'] = float(rt[1])
                    elif rt[0]=='End 2 Mismatch Rate':
                        alignqc_stats['R2_mismatch_rate'] = float(rt[1])
                    elif rt[0]=='High Quality Rate':
                        alignqc_stats['high_quality_rate'] = float(rt[1])
                    elif rt[0]=='Exonic Rate':
                        alignqc_stats['exon_mapping_rate'] = float(rt[1])
                    elif rt[0]=='Intronic Rate':
                        alignqc_stats['intron_mapping_rate'] = float(rt[1])
                    elif rt[0]=='Intergenic Rate':
                        alignqc_stats['intergenic_mapping_rate'] = float(rt[1])
                    elif rt[0]=='Intragenic Rate':
                        alignqc_stats['intragenic_mapping_rate'] = float(rt[1])
                    elif rt[0]=='Ambiguous Alignment Rate':
                        alignqc_stats['ambiguous_alignment_rate'] = float(rt[1])
                    elif rt[0]=='rRNA Rate':
                        alignqc_stats['rRNA_mapping_rate'] = float(rt[1])
                    elif rt[0]=='End 1 Sense Rate':
                        alignqc_stats['R1_sense_strand_mapping_rate'] = float(rt[1])
                    elif rt[0]=='End 2 Sense Rate':
                        alignqc_stats['R2_sense_strand_mapping_rate'] = float(rt[1])
                    elif rt[0]=='Avg. Splits per Read':
                        alignqc_stats['average_insertions_deletions_gaps_per_read'] = float(rt[1])
                    elif rt[0]=='Failed Vendor QC':
                        alignqc_stats['failed_vendor_QC'] = int(rt[1])
                    elif rt[0]=='Total Reads':
                        alignqc_stats['total_reads'] = int(rt[1])
                    elif rt[0]=='Low Mapping Quality':
                        alignqc_stats['low_mapping_quality'] = int(rt[1])
                    elif rt[0]=='Genes Detected':
                        alignqc_stats['genes_detected'] = int(rt[1])
                    elif rt[0]=="Mean 3' bias":
                        alignqc_stats['mean_gene_3p_bias_in_read_location'] = float(rt[1])
                    elif rt[0]=="Median 3' bias":
                        alignqc_stats['median_gene_3p_bias_in_read_location'] = float(rt[1])
                    elif rt[0]=="Median of Avg Transcript Coverage":
                        alignqc_stats['median_transcript_gene_coverage'] = float(rt[1])
    alignqc_stats['percent_failed_vendor_QC'] = 100.0*alignqc_stats['failed_vendor_QC']/alignqc_stats['total_reads'] if alignqc_stats['total_reads'] > 0 else 0
    alignqc_stats['percent_low_mapping_quality'] = 100.0*alignqc_stats['low_mapping_quality']/alignqc_stats['total_reads'] if	alignqc_stats['total_reads'] > 0 else 0
    
    return alignqc_stats


def summary_alignqc_rnaseqc( alignqc_stats ):
    """
    alignqc_stats = {'sample_file': '', 'mapping_rate_in_file': 0, 'unique_mapping_rate_of_mapped_in_file': 0, 'duplicate_mapping_rate_of_mapped_in_file': 0,
                     'base_mismatch_rate': 0, 'R1_mapping_rate': 0, 'R2_mapping_rate': 0, 'R1_mismatch_rate': 0, 'R2_mismatch_rate': 0,
                     'exon_mapping_rate': 0, 'intron_mapping_rate': 0, 'intergenic_mapping_rate': 0, 'intragenic_mapping_rate': 0,
                     'ambiguous_alignment_rate': 0, 'rRNA_mapping_rate': 0, 'high_qualiity_rate': 0,
                     'R1_sense_strand_mapping_rate': 0, 'R2_sense_strand_mapping_rate': 0,
                     'average_insertions_deletions_gaps_per_read': 0, 'failed_vendor_QC': 0, 'total_reads': 0, 'percent_failed_vendor_QC': 0,
                     'low_mapping_quality': 0, 'percent_low_mapping_quality': 0, 'genes_detected': 0,
                     'mean_gene_3p_bias_in_read_location': 0, 'median_gene_3p_bias_in_read_location': 0,
                     'median_transcript_gene_coverage': 0 }
    """
    alignqc_summary = {'percent_low_mapping_quality': {'result': 'PASS', 'suggest': '', 'value': 'NA', 'threshold': '<20%' },
                       'percent_failed_vendor_QC': {'result': 'PASS', 'suggest': '', 'value': 'NA', 'threshold': '<20%' },
                       'R1_mismatch_rate': {'result': 'PASS', 'suggest': '', 'value': 'NA', 'threshold': '<15%' },
                       'R2_mismatch_rate': {'result': 'PASS', 'suggest': '', 'value': 'NA', 'threshold': '<15%' },
                       'ambiguous_alignment_rate': {'result': 'PASS', 'suggest': '', 'value': 'NA', 'threshold': '<15%' },
                       'rRNA_mapping_rate': {'result': 'PASS', 'suggest': '', 'value': 'NA', 'threshold': '<20%' },
                       'average_insertions_deletions_gaps_per_read': {'result': 'PASS', 'suggest': '', 'value': 'NA', 'threshold': '<10%'}
                       }

    ## Report - Low mapping quality greater than 20%
    if float(alignqc_stats['percent_low_mapping_quality']) > 20.0:
        alignc_summary['percent_low_mapping_quality']['result'] = 'WARNING'
        alignc_summary['percent_low_mapping_quality']['suggest'] = 'Read Alignment - More than 20% of reads have low mapping quality. Consider removing these reads by setting an appropriate MAPQ threshold as a parameter during alignment.'
        alignc_summary['percent_low_mapping_quality']['suggest'] = alignqc_stats['percent_low_mapping_quality']
    
    ## Report - FAiled vendor QC > 20%
    if float(alignqc_stats['percent_failed_vendor_QC']) > 20.0:
        alignc_summary['percent_failed_vendor_QC']['result'] = 'FAIL'
        alignc_summary['percent_failed_vendor_QC']['suggest'] = 'Read Alignment - More than 20% of reads fail vendor QC. Consider re-sequencing this sample.'
        alignc_summary['percent_failed_vendor_QC']['value'] = alignqc_stats['percent_failed_vendor_QC']

    ## Report - R1 mismatch rate > 10%
    if float(alignqc_stats['R1_mismatch_rate']) > 10.0:
        alignc_summary['R1_mismatch_rate']['result'] = 'WARNING'
        alignc_summary['R1_mismatch_rate']['suggest'] = 'Read Alignment - R1 read mismatch rate > 10%. This should not affect gene expression results but may affect variant calling.'

    ## Report - R2 mismatch rate > 10%
    if float(alignqc_stats['R2_mismatch_rate']) > 10.0:
        alignc_summary['R2_mismatch_rate']['result'] = 'WARNING'
        alignc_summary['R2_mismatch_rate']['suggest'] = 'Read Alignment - R2 read mismatch rate > 10%. This should not affect gene expression results but may affect variant calling.'
        alignc_summary['R2_mismatch_rate']['value'] = alignqc_stats['R2_mismatch_rate']
        
    ## Report - ambiguous_alignment rate > 15%    
    if float(alignqc_stats['ambiguous_alignment_rate']) > 15.0:
        alignc_summary['ambiguous_alignment_rate']['result'] = 'WARNING'
        alignc_summary['ambiguous_alignment_rate']['suggest'] = 'Read Alignment - Ambiguous alignment rate > 15%. This may cause some error in gene read counts if ambiguous alignments are occuring in gene exon regions.'
        alignc_summary['ambiguous_alignment_rate']['value'] = alignqc_stats['ambiguous_alignment rate']
        
    ## Report - rRNA_mapping_rate > 20%
    if float(alignqc_stats['rRNA_mapping_rate']) > 20.0:
        alignc_summary['rRNA_mapping_rate']['result'] = 'WARNING'
        alignc_summary['rRNA_mapping_rate']['suggest'] = 'Read Alignment - Ribosomal RNA mapping rate > 20%. Consider ribosomal RNA removal in library prep.'
        alignc_summary['rRNA_mapping_rate']['value'] = alignqc_stats['rRNA_mapping_rate']
        
    ## Report - average_insertions_deletions_gaps_per_read > 10%
    if float(alignqc_stats['average_insertions_deletions_gaps_per_read']) > 10.0:
        alignc_summary['average_insertions_deletions_gaps_per_read']['result'] = 'WARNING'
        alignc_summary['average_insertions_deletions_gaps_per_read']['suggest'] = 'Read Alignment - Detected a large % of gaps, insertions and deletions (>10%).'
        alignc_summary['average_insertions_deletions_gaps_per_read']['value'] = alignqc_stats['average_insertions_deletions_gaps_per_read']
    
    return alignqc_summary


def getstats_fastqc( fastqc_files ):
    fastqc_stats_json = {}
    fastqc_summary_json = {}
    os.mkdir('../fastqc_temp')
    for sample_id, fastqc_file_list in fastqc_files.items():
        fastqc_stats_json[sample_id] = {}
        fastqc_summary_json[sample_id] = {}
        for j in range(0,len(fastqc_file_list)):
            fastqc_file = fastqc_file_list[j]
            with zipfile.ZipFile(os.path.join('../', fastqc_file),'r') as zfile:
                zfile.extractall('../fastqc_temp/')
            if os.path.exists('../fastqc_temp/{}/fastqc_data.txt'.format(fastqc_file[:-4])):
                if 'R2' in fastqc_file:
                    fastqc_stats_json[sample_id]['R2'] = parse_fastqc_fastqc('../fastqc_temp/{}/fastqc_data.txt'.format(fastqc_file[:-4]))
                    fastqc_summary_json[sample_id]['R2'] = summary_fastqc_fastqc(fastqc_stats_json[sample_id]['R2'])
                else:
                    fastqc_stats_json[sample_id]['R1'] = parse_fastqc_fastqc('../fastqc_temp/{}/fastqc_data.txt'.format(fastqc_file[:-4]))
                    fastqc_summary_json[sample_id]['R1'] = summary_fastqc_fastqc(fastqc_stats_json[sample_id]['R1'])
    return fastqc_stats_json, fastqc_summary_json


def getstats_aligner_rnastar( aligner_files ):
    aligner_stats_json, aligner_summary_json = {}, {}
    for sample_id, aligner_file in aligner_files.items():
        if os.path.exists( '../{}'.format(aligner_file)):
            aligner_stats_json[sample_id] = parse_aligner_rnastar('../{}'.format(aligner_file))
            aligner_summary_json[sample_id] = summary_aligner_rnastar(aligner_stats_json[sample_id])
            
    return aligner_stats_json, aligner_summary_json


def getstats_alignqc_rnaseqc( alignqc_files ):
    alignqc_stats_json, alignqc_summary_json = {}, {}
    for sample_id, alignqc_file in alignqc_files.items():
        if os.path.exists( '../{}'.format(alignqc_file)):
            alignqc_stats_json[sample_id] = parse_alignqc_rnaseqc('../{}'.format(alignqc_file))
            alignqc_summary_json[sample_id] = summary_alignqc_rnaseqc(alignqc_stats_json[sample_id])
    
    return alignqc_stats_json, alignqc_summary_json


def createSummaryTable_fastqc( fastqc_stats_json, outfile ):
    with open(outfile,'w') as fout:
        fout.write(','.join(['Sample_ID','Pair','Total_Sequences','Read_Length','Avg_Seq_Error_Rate','Peak_GC%'])+'\n')
        for sample_id, stats_full in fastqc_stats_json.items():
            for readpair, _stats in stats_full.items():
                total_raw_reads = str(_stats['total_raw_reads']) if 'total_raw_reads' in _stats else 'NA'
                read_length = str(_stats['read_length']) if 'read_length' in _stats else 'NA'
                avg_seq_error_rate = str(round(10**(-1.0*np.mean(_stats['per_base_quality']['mean'])/10),3))+'%' if 'per_base_quality' in _stats and 'mean' in _stats['per_base_quality'] else 'NA'
                max_gc_index = _stats['per_base_GC_content']['count'].index(max(_stats['per_base_GC_content']['count'])) if 'per_base_GC_content' in _stats and 'count' in _stats['per_base_GC_content'] else 'NA'
                max_gc = str(_stats['per_base_GC_content']['pos'][max_gc_index])+'%' if 'per_base_GC_content' in _stats and 'pos' in _stats['per_base_GC_content'] else 'NA'
                fout.write(','.join([sample_id, readpair, total_raw_reads, read_length, avg_seq_error_rate, max_gc])+'\n')
    return


def createSummaryTable_align( aligner_stats_json, alignqc_stats_json, outfile ):
    with open(outfile,'w') as fout:
        fout.write(','.join(['Input_Reads','%_Mapped','%_Uniquely_Mapped','%_Mismatch_Rate','%_Deletion_Rate','%_Insertion_Rate','%_rRNA_Reads', '%_Low_Mapping_Quality', '%_Failed_Vendor_QC', 'Genes_Detected'])+'\n')
        for sample_id in aligner_stats_json.keys():
            input_reads = str(aligner_stats_json[sample_id]['input_reads']) if sample_id in aligner_stats_json and 'input_reads' in aligner_stats_json[sample_id] else 'NA'
            p_mapped = str(aligner_stats_json[sample_id]['percent_mapped']) if sample_id in aligner_stats_json and 'percent_mapped' in aligner_stats_json[sample_id] else 'NA'
            p_unique_mapped = str(aligner_stats_json[sample_id]['percent_uniquely_mapped']) if sample_id in aligner_stats_json and 'percent_uniquely_mapped' in aligner_stats_json[sample_id] else 'NA'
            p_mm_rate = str(aligner_stats_json[sample_id]['percent_mismatch_rate']) if sample_id in aligner_stats_json and 'percent_mismatch_rate' in aligner_stats_json[sample_id] else 'NA'
            p_del_rate = str(aligner_stats_json[sample_id]['percent_deletion_rate']) if sample_id in aligner_stats_json and 'percent_deletion_rate' in aligner_stats_json[sample_id] else 'NA'
            p_ins_rate = str(aligner_stats_json[sample_id]['percent_insertion_rate']) if sample_id in aligner_stats_json and 'percent_insertion_rate' in aligner_stats_json[sample_id] else 'NA'
            p_rrna_rate = str(alignqc_stats_json[sample_id]['rRNA_mapping_rate']) if sample_id in alignqc_stats_json and 'rRNA_mapping_rate' in alignqc_stats_json[sample_id] else 'NA'
            p_low_mapq = str(alignqc_stats_json[sample_id]['low_mapping_quality']) if sample_id in alignqc_stats_json and 'low_mapping_quality' in alignqc_stats_json[sample_id] else 'NA'
            p_fail_qc = str(alignqc_stats_json[sample_id]['percent_failed_vendor_QC']) if sample_id in alignqc_stats_json and 'percent_failed_vendor_QC' in alignqc_stats_json[sample_id] else 'NA'
            p_genes = str(alignqc_stats_json[sample_id]['genes_detected']) if sample_id in alignqc_stats_json and 'genes_detected' in alignqc_stats_json[sample_id] else 'NA'

            fout.write(','.join( [input_reads,p_mapped,p_unique_mapped,p_mm_rate,p_del_rate,p_ins_rate,p_rrna_rate,p_low_mapq,p_fail_qc,p_genes] )+'\n')
    return


def addToSummary( summary_json, summary_cells, detail_cells ):
    """ add a summary JSON to summary reports
    """
    for sample, sstats in summary_json.items():
        for _stat, stat_result in sstats.items():
            # e.g., percent_mapped, {'result': 'PASS', 'suggest': . 'value':..} ..'sample'
            detail_cells[sample][_stat] = stat_result
            if _stat not in summary_cells:
                summary_cells[_stat] = stat_result
                summary_cells[_stat]['samples'] = [sample]
            elif stat_result['result'] == summary_cells[_stat]['result']:
                summary_cells[_stat]['samples'] += [sample]
            # highlight WARNING in summary report
            elif stat_result['result'] == 'WARNING' and summary_cells[_stat]['result'] == 'PASS':
                summary_cells[_stat] = stat_result
                summary_cells[_stat]['samples'] = [sample]
            elif stat_result['result'] == 'FAIL' and (summary_cells[_stat]['result'] in ['PASS', 'WARNING']):
                # highlight FAIL above all else
                summary_cells[_stat] = stat_result
                summary_cells[_stat]['samples'] = [sample]
    
    return summary_cells, detail_cells
            

def removeUnderscore( s ):
    """ Replace underscores with spaces, and capitalize
    """
    s_list = s.split('_')
    snew = ''
    for e in s_list:
        the_upper = e[0].upper() if 'RNA' not in e.upper() else e[0]
        snew += the_upper + e[1:] + ' '
    return snew


def createSummaryReports( fastqc_summary_json, aligner_summary_json, alignqc_summary_json, summary_html_file, details_html_file ):

    summary_cells = {} # per QC test
    detail_cells = {}  # per sample

    # get merge of all sample IDs
    for sample in fastqc_summary_json.keys():
        detail_cells[sample] = {}
    for sample in aligner_summary_json.keys():
        detail_cells[sample] = {}        
    for sample in alignqc_summary_json.keys():
        detail_cells[sample] = {}        

    summary_cells, detail_cells = addToSummary( fastqc_summary_json, summary_cells, detail_cells )
    summary_cells, detail_cells = addToSummary( aligner_summary_json, summary_cells, detail_cells )
    summary_cells, detail_cells = addToSummary( alignqc_summary_json, summary_cells, detail_cells )    

    summary_table = '<table border="1" style="padding:10px;text-align:left;"><tr><th>QC Stat</th><th>Status</th><th>Pass Criteria</th><th>Samples</th><th>Details</th></tr>'    
    for _stat, details in summary_cells.items():
        summary_table += '<tr>'
        summary_table += '<td>'+removeUnderscore(str(_stat))+'</td>'
        if details['result'] == 'PASS':
            result_style = '<td style="background-color:green;">'
            samples = 'ALL'
        elif details['result'] == 'WARNING':
            result_style = '<td style="background-color:yellow;">'
            samples = str(details['samples'])
        elif details['result'] == 'FAIL':
            result_style = '<td style="background-color:red;">'
            samples = str(details['samples'])
        else:
            result_style = '<td>'
            samples = str(details['samples'])
        summary_table += result_style+str(details['result'])+'</td>'
        summary_table += '<td>'+str(details['threshold'])+'</td>'
        summary_table += '<td>'+samples+'</td>'
        summary_table += '<td>'+str(details['suggest'])+'</td>'
        summary_table += '</tr>'
    summary_table += '</table>'
    summary_html = '<html><head><h1>RNA-Seq QC Summary</h1></head><body>'+summary_table+'</body></html>'

    details_tables = ''
    for sample, sstats in detail_cells.items():
        details_tables += '<h3>'+sample+'</h3><table border="1" style="padding:10px;text-align:left;"><tr><th>QC Stat</th><th>Status</th><th>Pass Criteria</th><th>Value</th><th>Details</th></tr>'
        for _stat, details in sstats.items():
            details_tables += '<tr>'
            details_tables += '<td>'+removeUnderscore(str(_stat))+'</td>'
            if details['result'] == 'PASS':
                result_style = '<td style="background-color:green;">'
            elif details['result'] == 'WARNING':
                result_style = '<td style="background-color:yellow;">'
            elif details['result'] == 'FAIL':
                result_style = '<td style="background-color:red;">'
            else:
                result_style = '<td>'            
            details_tables += result_style+str(details['result'])+'</td>'
            details_tables += '<td>'+str(details['threshold'])+'</td>'            
            details_tables += '<td>'+str(details['value'])+'</td>'            
            details_tables += '<td>'+str(details['suggest'])+'</td>'
            details_tables += '</tr>'
        details_tables += '</table>'
    details_html = '<html><head><h1>RNA-Seq QC details</h1></head><body>'+details_tables+'</body></html>'

    with open( summary_html_file, 'w' ) as fout:
        fout.write(summary_html)

    with open( details_html_file, 'w' ) as fout:
        fout.write(details_html)
    
    return


def run_program( arg_list ):
    """
    Parameters:
    -samples        comma-separated string of sample IDs to provide summary QC for
    -fastqc         fastq QC program used - default: fastqc
    -alignqc        alignment QC program used - default: rnaseqc
    -aligner        aligner program used - default: rnastar
    -input_fastqc   input fastq QC folder
    -input_alignqc  input align QC folder
    -input_aligner  input aligner folder
    -out            summary output folder
    """
    print('ARG LIST: {}'.format(str(arg_list)))
    input_fastqc = module_utils.getArgument( arg_list, '-input_fastqc', 'implicit' )  # input folders
    input_aligner = module_utils.getArgument( arg_list, '-input_aligner', 'implicit' )
    input_alignqc = module_utils.getArgument( arg_list, '-input_alignqc', 'implicit' )
    
    fastqc_program = module_utils.getArgument( arg_list, '-fastqc', 'implicit', 'fastqc' )  # programs used
    aligner_program = module_utils.getArgument( arg_list, '-aligner', 'implicit', 'rnastar' )
    alignqc_program = module_utils.getArgument( arg_list, '-alignqc', 'implicit', 'rnaseqc' )
    
    output_dir = module_utils.getArgument( arg_list, '-out' )  # local output dir
    sample_ids =  str(module_utils.getArgument( arg_list, '-samples', 'implicit', '' )).split(',')
    
    # get relevant files from programs - putting them all in one folder, hopefully file names don't clash
    os.chdir(output_dir)    
    all_input_files = os.listdir('../')
    fastqc_files = {}
    aligner_files = {}
    alignqc_files = {}
    # could probably do this in one loop, will optimize later
    for sample_id in sample_ids:
        for input_file in all_input_files:
            if sample_id in input_file and input_file.endswith('fastqc.zip') and str(fastqc_program).lower()=='fastqc':
                fastqc_files[sample_id] = fastqc_files[sample_id] + [input_file] if sample_id in fastqc_files else [input_file]
            elif sample_id in input_file and input_file.endswith('.Log.final.out') and (str(aligner_program).lower()=='rnastar' or str(aligner_program).lower()=='star'):
                aligner_files[sample_id] = input_file
            elif sample_id in input_file and input_file.endswith('.metrics.tsv') and str(alignqc_program).lower()=='rnaseqc':
                alignqc_files[sample_id] = input_file
    print('SAMPLES: '+str(sample_ids))
    print('FASTQC FILES: '+str(fastqc_files))
    print('ALIGNER FILES: '+str(aligner_files))
    print('ALIGNQC FILES: '+str(alignqc_files))    
    
    # now process each program, for each sample, to get summary stats
    fastqc_stats_json, fastqc_summary_json, aligner_stats_json, aligner_summary_json, alignqc_stats_json, alignqc_summary_json = {}, {}, {}, {}, {}, {}
    if fastqc_program=='fastqc':
        fastqc_stats_json, fastqc_summary_json = getstats_fastqc( fastqc_files )
        file_utils.writeJSON( fastqc_stats_json, 'fastqc_stats.json' )    
        file_utils.writeJSON( fastqc_summary_json, 'fastqc_summary.json' )
        createSummaryTable_fastqc( fastqc_stats_json, 'fastqc_summary.csv' )
        
    if aligner_program=='rnastar':
        aligner_stats_json, aligner_summary_json = getstats_aligner_rnastar( aligner_files )
        file_utils.writeJSON( aligner_stats_json, 'aligner_stats.json' )    
        file_utils.writeJSON( aligner_summary_json, 'aligner_summary.json' )
    
    if alignqc_program=='rnaseqc':
        alignqc_stats_json, alignqc_summary_json = getstats_alignqc_rnaseqc( alignqc_files )
        file_utils.writeJSON( alignqc_stats_json, 'alignqc_stats.json' )    
        file_utils.writeJSON( alignqc_summary_json, 'alignqc_summary.json' )
        
        createSummaryTable_align( aligner_stats_json, alignqc_stats_json, 'align_summary.csv' )

    createSummaryReports( fastqc_summary_json, aligner_summary_json, alignqc_summary_json, 'rnaseq_qc_report_summary.html', 'rnaseq_qc_report_sample_details.html' )
    
    return


if __name__ == '__main__':
    run_program( sys.argv[1:] )
