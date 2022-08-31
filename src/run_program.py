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
    fastqc_summary = {'per_base_quality': {'result': 'PASS', 'failed_pos_low': [], 'suggest': [] },
                      'per_base_content': {'result': 'PASS', 'failed_pos': [], 'failed_pos_low': [], 'failed_pos_high': [], 'failed_pos_low_nt': [], 'failed_pos_high_nt': [], 'suggest': [] },
                      'per_base_GC_content': {'result': 'PASS', 'suggest': [] },
                      'per_base_N_content': {'result': 'PASS', 'failed_pos': [], 'suggest': [] },
                      'deduplicated_perc': {'result': 'PASS', 'suggest': [] },
                      'adapter_content': {'result': 'PASS', 'failed_pos': [], 'suggest': [] }
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
            fastqc_summary['per_base_quality']['suggest'].append("Read Sequencing - low sequencing quality detected on the 5'end. Suggest trimming {} bases off the 5' end of every read.".format(str(current_pos)))

        pos_len = len(fastqc_summary['per_base_quality']['failed_pos_low'])
        max_pos = max(fastqc_stats['per_base_quality']['pos'])
        current_pos = max_pos
        for i in range(0,pos_len):
            if fastqc_summary['per_base_quality']['failed_pos_low'][pos_len-i-1] == current_pos:
                current_pos = current_pos - 1
            else:
                break
        if current_pos < max_pos:
            fastqc_summary['per_base_quality']['suggest'].append("Read Sequencing - low sequencing quality detected on the 5'end. Suggest trimming {} bases off the 3' end of every read.".format(str(max_pos-current_pos)))

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
            fastqc_summary['per_base_content']['suggest'].append("Read sequencing - detected nucleotide bias on the 5'end. Suggest trimming {} bases off the 5' end of every read.".format(str(current_pos)))

        pos_len = len(fastqc_summary['per_base_content']['failed_pos_low'])
        max_pos = max(fastqc_stats['per_base_content']['pos'])
        current_pos = max_pos
        for i in range(0,pos_len):
            if fastqc_summary['per_base_content']['failed_pos_low'][pos_len-i-1] == current_pos:
                current_pos -=1
            else:
                break
        if current_pos < max_pos:
            fastqc_summary['per_base_content']['suggest'].append("Read sequencing - detected nucleotide bias on the 3'end. Suggest trimming {} bases off the 3' end of every read.".format(str(max_pos-current_pos)))

    ## report extreme GC% bias
    max_GC_perc_index = fastqc_stats['per_base_GC_content']['count'].index(max(fastqc_stats['per_base_GC_content']['count']))
    max_GC_perc = fastqc_stats['per_base_GC_content']['pos'][max_GC_perc_index]
    if max_GC_perc < 20 or max_GC_perc > 80:
        fastqc_summary['per_base_GC_content']['result'] = 'FAIL'
        fastqc_summary['per_base_GC_content']['suggest'] = 'Read sequencing - detected extreme GC bias - peak GC% is {}%'.format(str(max_GC_perc))
    elif max_GC_perc < 30 or max_GC_perc > 70:
        fastqc_summary['per_base_GC_content']['result'] = 'WARNING'
        fastqc_summary['per_base_GC_content']['suggest'] = 'Read sequencing - detected GC bias - peak GC% is {}%'.format(str(max_GC_perc))

    ## report any positions that have a high number of Ns
    max_pos = max(fastqc_stats['per_base_N_content']['pos'])    
    for i in range(0,len(fastqc_stats['per_base_N_content']['pos'])):
        if float(fastqc_stats['per_base_N_content']['perc'][i]) > 10:
            fastqc_summary['per_base_N_content']['failed_pos'].append(int(fastqc_stats['per_base_N_content']['pos'][i]))
    suggested_mismatch_rate = int(100.0*len(fastqc_summary['per_base_N_content']['failed_pos'])/max_pos)
    if len(fastqc_summary['per_base_N_content']['failed_pos']) > 0 and suggested_mismatch_rate <= 15:
        fastqc_summary['per_base_N_content']['result'] = 'WARNING'
        fastqc_summary['per_base_N_content']['suggest'] = 'Read sequencing - detected greater than 10% unknown base calls (N) at positions {}. Suggest increasing allowed mismatch rate to {}% during alignment.'.format(len(fastqc_summary['per_base_N_content']['failed_pos']), suggested_mismatch_rate )
    elif len(fastqc_summary['per_base_N_content']['failed_pos']) > 0 and suggested_mismatch_rate > 15:
        fastqc_summary['per_base_N_content']['result'] = 'FAIL'
        fastqc_summary['per_base_N_content']['suggest'] = 'Read sequencing - detected greater than 10% unknown base calls (N) at positions {}. Too many N positions - alignment would require an allowed mismatch rate of greater than 15%.'.format(str(len(fastqc_summary['per_base_N_content']['failed_pos'])), str(suggested_mismatch_rate) )

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
            fastqc_summary['adapter_content']['suggest'].append("Read sequencing - detected sequencing adapter on the 5'end. Suggest trimming {} bases off the 5' end of every read.".format(str(current_pos)))

        pos_len = len(fastqc_summary['adapter_content']['failed_pos'])
        max_pos = max(fastqc_stats['adapter_content']['pos'])
        current_pos = max_pos
        for i in range(0,pos_len):
            if fastqc_summary['adapter_content']['failed_pos_low'][pos_len-i-1] == current_pos:
                current_pos -=1
            else:
                break
        if current_pos < max_pos:
            fastqc_summary['adapter_content']['suggest'].append("Read sequencing - detected sequencing adapter on the 3'end. Suggest trimming {} bases off the 3' end of every read.".format(str(max_pos-current_pos)))
    
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
    aligner_summary = {'percent_mapped': {'result': 'PASS', 'suggest': [] },
                       'percent_uniquely_mapped': {'result': 'PASS', 'suggest': [] },
                       'percent_mismatch_rate': {'result': 'PASS', 'suggest': [] },
                       'percent_deletion_rate': {'result': 'PASS', 'suggest': [] },
                       'percent_insertion_rate': {'result': 'PASS', 'suggest': [] }
                      }

    print('ALIGNER STATS: '+str(aligner_stats))
    
    ## report if % mapped < 50%
    if float(aligner_stats['percent_mapped']) < 50.0:
        aligner_summary['percent_mapped']['result'] = 'WARNING'
        aligner_summary['percent_mapped']['suggest'] = 'Alignment mapping rate is low [{}%]. Consider read trimming or increasing allowed mismatch rate as suggested through FASTQ QC, or ignoring this sample if possible.'.format( str(round(float(aligner_stats['percent_mapped']),2)))
    elif float(aligner_stats['percent_mapped']) < 20.0:
        aligner_summary['percent_mapped']['result'] = 'FAIL'
        aligner_summary['percent_mapped']['suggest'] = 'Alignment mapping rate is low [{}%]. Consider removing this sample from further analysis.'.format( str(round(float(aligner_stats['percent_mapped']),2)))

    ## report if % uniquely mapped is < 50% of mapped
    if float(aligner_stats['percent_mapped']) > 0 and float(aligner_stats['percent_uniquely_mapped'])/float(aligner_stats['percent_mapped']) < 0.5:
        aligner_summary['percent_uniquely_mapped']['result'] = 'WARNING'
        aligner_summary['percent_uniquely_mapped']['suggest'] = 'More than 50% of mapped reads are aligning to multiple regions. If this is an issue, consider more stringent mapping parameters.'

    ## report if mismatch rate > 15%
    if float(aligner_stats['percent_mismatch_rate']) > 15.0:
        aligner_summary['percent_mismatch_rate']['result'] = 'WARNING'
        aligner_summary['percent_mismatch_rate']['suggest'] = 'Mismatch rate is {}%. Consider more stringent mapping parameters, or re-sequencing if possible.'.format(str(round(float(aligner_stats['percent_mismatch_rate']),2)))

    ## report if deletion rate > 5%
    if float(aligner_stats['percent_deletion_rate']) > 5.0:
        aligner_summary['percent_deletion_rate']['result'] = 'WARNING'
        aligner_summary['percent_deletion_rate']['suggest'] = 'Deletion rate is {}%. Consider more stringent mapping parameters, or re-sequencing if possible.'.format(str(round(float(aligner_stats['percent_deletion_rate']),2)))

    ## report if insertion rate > 5%
    if float(aligner_stats['percent_insertion_rate']) > 5.0:
        aligner_summary['percent_insertion_rate']['result'] = 'WARNING'
        aligner_summary['percent_insertion_rate']['suggest'] = 'Insertion rate is {}%. Consider more stringent mapping parameters, or re-sequencing if possible.'.format(str(round(float(aligner_stats['percent_insertion_rate']),2)))        
    
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
    alignqc_summary = {'percent_low_mapping_quality': {'result': 'PASS', 'suggest': '' },
                       'percent_failed_vendor_QC': {'result': 'PASS', 'suggest': '' },
                       'R1_mismatch_rate': {'result': 'PASS', 'suggest': '' },
                       'R2_mismatch_rate': {'result': 'PASS', 'suggest': '' },                       
                       'ambiguous_alignment_rate': {'result': 'PASS', 'suggest': '' },
                       'rRNA_mapping_rate': {'result': 'PASS', 'suggest': '' },
                       'average_insertions_deletions_gaps_per_read': {'result': 'PASS', 'suggest': '' }                       
                       }

    ## Report - Low mapping quality greater than 20%
    if float(alignqc_stats['percent_low_mapping_quality']) > 20.0:
        alignc_summary['percent_low_mapping_quality']['result'] = 'WARNING'
        alignc_summary['percent_low_mapping_quality']['suggest'] = 'Read Alignment - {}% of reads are of low mapping quality. Consider removing these reads by setting an appropriate MAPQ threshold as a parameter during alignment.'.format(str(alignqc_stats['percent_low_mapping_quality']))

    ## Report - FAiled vendor QC > 20%
    if float(alignqc_stats['percent_failed_vendor_QC']) > 20.0:
        alignc_summary['percent_failed_vendor_QC']['result'] = 'FAIL'
        alignc_summary['percent_failed_vendor_QC']['suggest'] = 'Read Alignment - {}% of reads fail vendor QC. Consider re-sequencing this sample.'.format(str(alignqc_stats['percent_failed_vendor_QC']))

    ## Report - R1 or R2 mismatch rate > 10%
    if float(alignqc_stats['R1_mismatch_rate']) > 10.0:
        alignc_summary['R1_mismatch_rate']['result'] = 'WARNING'
        alignc_summary['R1_mismatch_rate']['suggest'] = 'Read Alignment - R1 read mismatch rate is {}%. This should not affect gene expression results but may affect variant calling.'

    ## Report - ambiguous_alignment rate > 15%    
    if float(alignqc_stats['R2_mismatch_rate']) > 15.0:
        alignc_summary['R2_mismatch_rate']['result'] = 'WARNING'
        alignc_summary['R2_mismatch_rate']['suggest'] = 'Read Alignment - Ambiguous alignment rate is {}%. This may cause some error in gene read counts if ambiguous alignments are occuring in gene exon regions.'

    ## Report - rRNA_mapping_rate > 20%
    if float(alignqc_stats['rRNA_mapping_rate']) > 20.0:
        alignc_summary['rRNA_mapping_rate']['result'] = 'WARNING'
        alignc_summary['rRNA_mapping_rate']['suggest'] = 'Read Alignment - Ribosomal RNA mapping rate is {}%. Consider ribosomal RNA removal in library prep.'.format(str(alignqc_stats['rRNA_mapping_rate']))
        
    ## Report - average_insertions_deletions_gaps_per_read > 10%
    if float(alignqc_stats['average_insertions_deletions_gaps_per_read']) > 15.0:
        alignc_summary['average_insertions_deletions_gaps_per_read']['result'] = 'WARNING'
        alignc_summary['average_insertions_deletions_gaps_per_read']['suggest'] = 'Read Alignment - Detected a large % of gaps, insertions and deletions (total {}%).'.format(str(alignqc_stats['average_insertions_deletions_gaps_per_read']))
    
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
    
    return


if __name__ == '__main__':
    run_program( sys.argv[1:] )
