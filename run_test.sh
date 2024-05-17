set -ex

cd test/input/fastq/

# bash pull-data.sh

cd -

nextflow run main.nf -profile test

bash md5.sh
