aws s3 cp $1.test.job.json s3://hubseq-data/modules/$1/job/
aws s3 cp $1.test.io.json s3://hubseq-data/modules/$1/io/
aws s3 cp $1.dryrun_test.io.json s3://hubseq-data/modules/$1/io/
aws s3 cp $1.template.json s3://hubseq-data/templates/
