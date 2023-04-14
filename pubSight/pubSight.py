# importing libraries

import pandas as pd
from .utils import *
import argparse
import warnings
import os


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', help="path to file that has all the query terms",
                        type=str, required=False)
    parser.add_argument('--email', '-e', help="a valid email address",
                        type=str, required=True)
    parser.add_argument('--out_dir', '-o', help="path to report directory",
                        type=str, required=True)
    parser.add_argument('--color_palette', '-c', help="matplotlib color palette",
                        type=str, required=False, default='cividis')
    parser.add_argument('--group_legend', '-g', help="make a group legend for whole plot",
                        action='store_true', default=False)
    parser.add_argument('--pubmed_data', '-p', help="path to a csv dataframe that has four columns"
                                                    " ['year', 'count', 'main_term', 'sub_term']", type=str,
                        required=False)

    return parser.parse_args()
def check_requirements(args):
    """
    Check requirements (file format, dependencies, permissions)
    """
    # check the third party softwares for plotting the results
    try:
        import pandas as pd
    except ImportError:
        sys.exit("--- Please check your installation for pandas library")
    # Check that the output directory is writeable
    output_dir = os.path.abspath(args.out_dir)
    if not os.path.isdir(output_dir):
        try:
            print("Creating output directory: " + output_dir)
            os.mkdir(output_dir)
        except EnvironmentError:
            sys.exit("CRITICAL ERROR: Unable to create output directory.")



def main():
    # Parse arguments from command line
    warnings.filterwarnings("ignore")
    warnings.simplefilter('ignore')

    args = parse_arguments()

    print(args)  # printing Namespace
    check_requirements(args)
    assert (args.input is not None) or (args.pubmed_data is not None), "You should provide input queries or pubmed data"
    if args.pubmed_data is not None:
        print('Loading data from local memory')
        df_main = read_terms(args.pubmed_data, delimiter='\t')
        print(df_main)
    else:
        print('Fetching data from pubmed')
        df = read_terms(args.input, delimiter='\t')
        df_main = get_from_pd(data=df, year=2023, email=args.email, write=True, report_dir=args.out_dir)

    pubmed_plot(data=df_main, colormap=args.color_palette, group_legend=args.group_legend, report_dir=args.out_dir)
    return print('done!')


# main()
if __name__ == "__main__":
    main()
