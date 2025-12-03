import argparse
import pandas as pd

def main(**kwargs):

    df = pd.read_csv(kwargs['path'], sep="\t")
    df = df.drop_duplicates(subset=['alteration_type','alteration','treatment','tumorType'],ignore_index=True)
    fn = kwargs['path'].split('.')[0]
    df.to_csv(fn+"_dedup.csv", sep="\t")

if __name__ == "__main__":

    def add_arguments(parser):
        parser.add_argument('--path', type=str, help='')

    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()
    main(**vars(args))
