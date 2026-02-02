
import argparse
import pandas as pd


def make_clade_file(input_file, output_file, clade):
    # read the master qc file
    qc_df = pd.read_csv(input_file)
    # subset to only the matching clade
    clade_df = qc_df[qc_df['auriclass_clade'] == clade]
    # write out a one-column file with only the isolate names
    clade_df[['Sample']].to_csv(output_file, index=False, header=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a master qc file.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide a file name for the output file.''',
        default=None
        )
    parser.add_argument(
        '--clade','-c',choices=['Clade I','Clade III','Clade IV'],
        help='''Specify a clade'''
        )    
    args = parser.parse_args()
    make_clade_file(args.input, args.output, args.clade)

if __name__ == '__main__':
    main()
