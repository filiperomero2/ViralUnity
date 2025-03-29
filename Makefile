.PHONY: test run-meta run-consensus install build-docker run-docker

test:
	python3 -m unittest discover ./test -p *test.py

install:
	pip install -e .

run-meta:
	viralunity meta \
		--data-type illumina \
		--sample-sheet input/viralunity_samplesheet.csv \
		--config-file output/config_meta.yml \
		--run-name test-meta \
		--kraken2-database input/database/kraken2 \
		--krona-database input/database/krona/taxonomy/ \
		--adapters input/references/SARS-CoV-2_RefSeq.fasta \
		--threads 6 \
		--threads-total 6 \
		--output output/test-meta-exmaple

run-consensus:
	viralunity consensus \
		--data-type illumina \
		--sample-sheet input/viralunity_samplesheet.csv \
		--config-file output/config_meta.yml \
		--run-name test-meta \
		--kraken2-database input/database/kraken2 \
		--krona-database input/database/krona/taxonomy/ \
		--adapters input/references/SARS-CoV-2_RefSeq.fasta \
		--threads 1 \
		--threads-total 1 \
		--output output/test-meta-exmaple

build-docker:
	docker build -t viralunity/viralunity:latest .

run-docker:
	docker run --rm -i -t viralunity/viralunity:latest
