import argparse
import sys
from viralunity import __program__, __version__
from viralunity.viralunity_meta import fill_arg_parser_meta
from viralunity.viralunity_consensus import fill_arg_parser_consensus

def cli(args):
    parser = argparse.ArgumentParser(
        description="ViralUnity is a simple tool to perform analysis of viral high-throughput sequencing data.",
        usage="""viralunity <tool> [options]""",
    )

    parser.add_argument('-v', '--version', action='version', version=f'{__program__} {__version__}')
    
    subparsers = parser.add_subparsers()
    meta_parser = subparsers.add_parser('meta', description="A script to generate config files and run the viralunity metagenomics pipeline")
    consensus_parser = subparsers.add_parser('consensus',  description="A script to generate config files and run the viralunity consensus pipeline")
    
    fill_arg_parser_meta(meta_parser)
    fill_arg_parser_consensus(consensus_parser)
    
    arguments = parser.parse_args(args)
    
    match arguments.tool:
        case 'meta':
            from viralunity.viralunity_meta import main as meta_main
            return meta_main(args)
        case 'consensus':
            from viralunity.viralunity_consensus import main as consensus_main
            return consensus_main(args)
        case _:
            return 1
    

def main():
    sys.exit(cli(sys.argv[1:]))