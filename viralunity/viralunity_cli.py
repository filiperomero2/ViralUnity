import argparse
import sys
from viralunity import __program__, __version__
from viralunity.viralunity_consensus_cli import fill_arg_parser_consensus
from viralunity.viralunity_meta import main as meta_main
from viralunity.viralunity_consensus import main as consensus_main
from viralunity.viralunity_meta_cli import fill_arg_parser_meta
    
def cli(args):
    viralunity_parser = argparse.ArgumentParser(
        description="ViralUnity is a simple tool to perform analysis of viral high-throughput sequencing data.",
    )

    viralunity_parser.add_argument('-v', '--version', action='version', version=f'{__program__} {__version__}')
    
    subparsers = viralunity_parser.add_subparsers()
    
    fill_arg_parser_meta(subparsers)
    fill_arg_parser_consensus(subparsers)
    
    arguments = viralunity_parser.parse_args(args)
    try:
        arguments.func(vars(arguments))
    except Exception as e:
        print(f"Error: {e}")
        viralunity_parser.print_help()
        raise e
        

def main():
    sys.exit(cli(sys.argv[1:]))

if __name__ == "__main__":
    main()