import os, subprocess, sys
import pandas as pd
import numpy as np
sys.path.append('global_utils/src/')
# sys.path.append('/global_utils/src/')
import module_utils
import aws_s3_utils
import file_utils

def run_program( arg_list ):
    """
    Parameters:
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
    
    # get relevant files from programs - hopefully file names don't clash
    print(os.listdir(os.getcwd()))
    print(' NOW TO OUTPUT DIR: ')
    os.chdir(output_dir)
    print(os.listdir('../'))
    return


if __name__ == '__main__':                                                                                                          
    run_program( sys.argv[1:] )
