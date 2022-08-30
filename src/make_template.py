import json
MODULE = 'rnaseq_summaryqc'

mi_template_json = {'module_version': '00.00.00', 'program_name': 'rnaseq_summaryqc', 'program_subname': '', 'program_version': '1.0', 'compute': {'environment': 'aws', 'language': 'Python', 'language_version': '3.7', 'vcpus': 2, 'memory': 8000}, 'program_arguments': '-fastqc fastqc -aligner rnastar -alignqc rnaseqc', 'program_input': [{'input_type': 'folder', 'input_file_type': '', 'input_position': 0, 'input_prefix': '-input_fastqc,-input_aligner,-input_alignqc'}], 'program_output': [{'output_type': 'folder', 'output_file_type': '', 'output_position': -1, 'output_prefix': '-out'}], 'alternate_inputs': [], 'alternate_outputs': [], 'defaults': {}}
with open(MODULE+'.template.json','w') as fout:
    json.dump(mi_template_json, fout)

io_dryrun_json = {'input': ['s3://hubtenants/test/rnaseq/run_test1/fastqc/,s3://hubtenants/test/rnaseq/run_test1/rnastar/,s3://hubtenants/test/rnaseq/run_test1/rnaseqc/'], 'output': ['s3://hubtenants/test/rnaseq/run_test1/rnaseq_summaryqc/'], 'alternate_inputs': [], 'alternate_outputs': [], 'program_arguments': '', 'sample_id': MODULE+'_test', 'dryrun': ''}
io_json = {'input': ['s3://hubtenants/test/rnaseq/run_test1/fastqc/,s3://hubtenants/test/rnaseq/run_test1/rnastar/,s3://hubtenants/test/rnaseq/run_test1/rnaseqc/'], 'output': ['s3://hubtenants/test/rnaseq/run_test1/rnaseq_summaryqc/'], 'alternate_inputs': [], 'alternate_outputs': [], 'program_arguments': '', 'sample_id': MODULE+'_test'}

with open(MODULE+'.dryrun_test.io.json','w') as fout:
    json.dump(io_dryrun_json, fout)
with open(MODULE+'.test.io.json','w') as fout:
    json.dump(io_json, fout)

# job info test JSONs                                                                                                        
job_json = {"container_overrides": {"command": ["--module_name", MODULE, "--run_arguments", "s3://hubseq-data/test/modules/"+MODULE+"/job/"+MODULE+".test.io.json", "--working_dir", "/home/"]}, "jobqueue": "batch_scratch_queue", "jobname": "job_"+MODULE+"_test"}
with open(MODULE+'.test.job.json','w') as fout:
    json.dump(job_json, fout)
