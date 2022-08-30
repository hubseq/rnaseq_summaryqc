#
# run_main
#
# Template wrapper script that runs a command-line program within a Docker container.
#
import os, subprocess, sys
from datetime import datetime
sys.path.append('global_utils/src/')
import module_utils

def runOtherPre( input_dir, output_dir, run_json ):
    """ This function is used to run any other commands BEFORE the main program has run.
    run_json has most of what you might need to run other commands, and has the following structure:

    run_json = {'module': module_name, 'local_input_dir': <LOCAL_INPUT_DIR>, 'local_output_dir': <LOCAL_OUT_DIR>, \
		'remote_input_dir': remote_input_directory, 'remote_output_dir': remote_output_directory, \
                'program_arguments': program_arguments, 'run_arguments': run_arguments_json, \
                'module_instance_json': module_instance_json}

    LOCAL_INPUT_DIR has any downloaded files. LOCAL_OUT_DIR has any output data or log files that will be uploaded.

    If you are not running any other commands or post-processing, then leave this function blank.
    """
    remote_output_dir = run_json['remote_output_dir']
    pargs = run_json['program_arguments']
    pargs_list = pargs.split(' ')
    pargs_list = ['python','/run_program.py','-remoteout',remote_output_dir] + pargs_list[1:]
    run_json['program_arguments'] = ' '.join(pargs_list)
    print('PRE LIST DIRS CWD: {}'.format(str(os.listdir(os.getcwd()))))
    print('PRE LIST DIRS ROOT: {}'.format(str(os.listdir('/'))))
    return run_json


def runOtherPost( input_dir, output_dir, run_json ):
    """ This function is used to run any other commands AFTER the main program has run.
    run_json has most of what you might need to run other commands, and has the structure shown above.

    If you are not running any other commands or pre-processing, then leave this function blank.
    """
    return run_json


def runMain():
    # time the duration of module container run
    run_start = datetime.now()
    print('Container running...')
    
    # initialize program run
    run_json = module_utils.initProgram()
    
    # do any pre-processing (specific to module)
    run_json = runOtherPre( run_json['local_input_dir'], run_json['local_output_dir'], run_json )
    
    # run main program
    print('RUN_JSON: {}'.format(str(run_json)))
    module_utils.runProgram( run_json['program_arguments'], run_json['local_output_file'] )
    
    # do any post-processing
    run_json = runOtherPost( run_json['local_input_dir'], run_json['local_output_dir'], run_json )

    # create run log that includes program run duration
    run_end = datetime.now()
    run_json['module_run_duration'] = str(run_end - run_start)
    module_utils.logRun( run_json, run_json['local_output_dir'] )
    
    # upload output data files
    module_utils.uploadOutput( run_json['local_output_dir'], run_json['remote_output_dir'] )
    print('DONE!')
    
    return


if __name__ == '__main__':
    runMain()
