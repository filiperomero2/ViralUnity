.PHONY: test run

test:
	python3 -m unittest discover ./test -p *test.py

run:
	viralunity_meta \
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