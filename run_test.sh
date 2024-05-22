set -ex

cd test/input/fastq/

# bash pull-data.sh

cd -

nextflow run main.nf -profile test,mamba

bash md5.sh
