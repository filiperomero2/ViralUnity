from setuptools import setup, find_packages
from viralunity import __version__, __program__, _description

setup(
    name='viralunity',
    version=__version__,
    packages=find_packages(),
    scripts=[
        'viralunity/scripts/consensus_illumina.smk',
        'viralunity/scripts/consensus_nanopore.smk',
        'viralunity/scripts/metagenomics_illumina.smk',
        'viralunity/scripts/metagenomics_nanopore.smk',
        'viralunity/scripts/metagenomics.smk',
        'viralunity/scripts/qc.smk',
    ],
    include_package_data=True,
    install_requires=[
        'pandas>=2.0.3',
        'biopython>=1.81',
        'snakemake',
        # Add other dependencies here
    ],
    description=_description,
    url='https://github.com/filiperomero2/ViralUnity',
    author='Filipe Moreira, Felippe Nacif',
    author_email='filiperomero2@gmail.com',
    entry_points={
        'console_scripts': [
            f'{__program__}_meta=viralunity.viralunity_meta:main',
            f'{__program__}_consensus=viralunity.viralunity_consensus:main',
            f'{__program__}_create_samplesheet=viralunity.viralunity_create_samplesheet:main',
        ],
    },
    keywords=[], # TODO: Add keywords @Filipe
    zip_safe=False,
)