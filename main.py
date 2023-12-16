import argparse
from src.preanalysis import Preanalysis


def my_flags():
    """
    Parse command line arguments using argparse.

    Returns:
        argparse.Namespace: An object containing the parsed arguments.
    """
    parser = argparse.ArgumentParser(description='This is the parser for the CCS project')

    # Add the -P or --preanalysis flag
    parser.add_argument('-P', '--preanalysis', action='store_true',
                        help='perform the preanalysis of the dataset')

    # Add the -N or --normalization flag
    parser.add_argument('-N', '--normalization', action='store_true',
                        help='use normalization during the preanalysis')

    return parser.parse_args()


if __name__ == '__main__':
    """
    Main entry point of the script.

    Parses command line arguments, creates an instance of the Preanalysis class,
    and performs pre-analysis based on the provided flags.
    """
    my_preanalysis = Preanalysis(my_flags().normalization) if my_flags().preanalysis else None
