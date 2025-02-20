# Description: Makefile for Nextflow pipeline development

.PHONY: clean md5sum test run-test

clean:
	rm -rf work/ testing/output .nextflow.log* .nextflow/

md5sum:
	bash assets/test/md5.sh

run-test:
	nextflow run main.nf -profile test,mamba

test: clean
	$(MAKE) run-test
	$(MAKE) md5sum
	# check if this file is not present testing/output/demixed/SRR30789620.demixed
	@if [ ! -f testing/output/demixed/SRR30789620.demixed ]; then \
		echo "File testing/output/demixed/SRR30789620.demixed is not present"; \
	else \
		echo "Error: File testing/output/demixed/SRR30789620.demixed should not exist"; \
		exit 1; \
	fi
