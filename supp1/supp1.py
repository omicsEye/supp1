# importing libraries

import pandas as pd
from supp1.utils import *
import argparse
import warnings


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--query_terms', '-qt', help="path to file that has all the query terms",
                        type=str, required=True)
    parser.add_argument('--email', '-em', help="a valid email address",
                        type=str, required=True)
    parser.add_argument('--out_dir', '-od', help="path to report directory",
                        type=str, required=True)
    parser.add_argument('--color_palette', '-cp', help="matplotlib color palette",
                        type=str, required=False, default='cividis')
    parser.add_argument('--pubmed_data', '-pd', help="path to a csv dataframe that has four columns"
                                                     " ['year', 'count', 'main_term', 'sub_term']", type=str,
                        required=False)

    return parser.parse_args()


def main():
    # Parse arguments from command line
    warnings.filterwarnings("ignore")
    warnings.simplefilter('ignore')

    args = parse_arguments()

    print(args)  # printing Namespace
    print('Fetching data from pubmed')
    if args.pubmed_data is not None:
        df_main = read_terms(args.pubmed_data, delimiter=',')
    else:
        df = read_terms(args.query_terms)
        df_main = get_from_pd(data=df, year=2023, email=args.email, write=True, report_dir=args.out_dir)

    pubmed_plot(data=df_main, colormap=args.color_palette, report_dir=args.out_dir)
    return print('done!')


# main()
if __name__ == "__main__":
    main()
