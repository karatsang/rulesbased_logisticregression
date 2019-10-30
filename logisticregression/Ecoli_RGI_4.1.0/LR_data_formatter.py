#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
import pandas as pd
import itertools
import numpy as np

def get_parser():
    """
    Create the argument parser
    """
    parser = argparse.ArgumentParser(description='Encode ')
    parser.add_argument('rgi_out', type=str, help='Folder containing RGI output')
    parser.add_argument('ast_tsv', type=str, help='Filepath to AST TSV')
    parser.add_argument('--perfect', dest='perfect_only', default=False,
                        action='store_true', help='Only use perfect RGI results')
    return parser


def get_rgi_df(rgi_folder):
    """
    Parse folder of RGI tsvs into a single large dataframe
    """
    rgi_results = []
    genome_ids = []
    for rgi_tsv in glob.glob(os.path.join(rgi_folder, '*.txt')):
        df = pd.read_csv(rgi_tsv, sep='\t')
        genome_id = os.path.basename(rgi_tsv.rstrip('.txt'))
        df['Sample'] = genome_id
        genome_ids.append(genome_id)
        rgi_results.append(df)

    # combine and clean up the rgi tsvs into one big tsv
    rgi_results = pd.concat(rgi_results)
    rgi_results = rgi_results.reset_index()

    return rgi_results, genome_ids


def encode_rgi(rgi_df, genome_ids):
    """
    Encode the RGI results into a matrix where each row is one 'genome' i.e.
    single RGI tsv ID and each column is one of the set of AROs present in
    any of the RGI tsvs in the input folder
    """
    rgi_encoded = pd.DataFrame(index=genome_ids,
                               columns=rgi_df['Best_Hit_ARO'].unique()).fillna(0)
    # print(rgi_encoded)
    for genome_id, rgi_data in rgi_df.iterrows():
        rgi_encoded.loc[rgi_data['Sample'], rgi_data['Best_Hit_ARO']] += 1

    return rgi_encoded


def gather_term_into_dict(rgi_df, genome_ids, column):
    """
    Parse the RGI df and based on column of choice e.g. AMR Gene Family,
    Drug Class or Mechanism return a dictionary mapping each genome ID
    to a list of resistances from that column.

    I.e. {genome_X: [all drug class resistances]}
    """
    data_dict = {x: [] for x in genome_ids}
    for row in rgi_df.fillna('').iterrows():
        data_dict[row[1]['Sample']] += row[1][column].split(';')
    return data_dict


def get_tag_encoded_df(rgi_df, genome_ids, column):
    """
    Get a ARO tag e.g. drug class and AMR gene family encoded dataframes from
    the combined RGI output dataframe

    Needs a separate version to normal encoding due to the list nature of
    some of the tags (i.e. multiple semi-colon separated for a single RGI
    result row)
    """
    class_dict = gather_term_into_dict(rgi_df, genome_ids, column)

    # get all values
    all_class_values = set(itertools.chain.from_iterable(class_dict.values()))

    class_df = pd.DataFrame(index=genome_ids, columns=all_class_values)
    class_df = class_df.fillna(0)
    for genome_id, class_list in class_dict.items():
        for class_value in class_list:
            class_df.loc[genome_id, class_value] += 1

    return class_df


if __name__ == '__main__':

    parser = get_parser()
    args = parser.parse_args()

    # read the AST data (assuming its formatted as a TSV with a column ID
    # which corresponds to ID from the RGI output TSV names)
    ast_df = pd.read_csv(args.ast_tsv, sep='\t')
    ast_df = ast_df.where(ast_df.notnull(), None)

    # encode lable R = 1, S = 0
    ast_df[ast_df == 'R'] = 1
    ast_df[ast_df == 'RA'] = 1
    ast_df[ast_df == 'I'] = 1   
    ast_df[ast_df == 'S'] = 0
    # ast_df = ast_df.where(ast_df.notnull(), None)
    ast_df = ast_df.set_index("Sample")
    # ast_df = ast_df.astype('int64')
    print("Saving AST to ast_df.pkl")
    ast_df.to_pickle('ast_df.pkl')

    # read the folder of RGI output TSVs into a single dataframe with an
    # ID determined from the name of the TSV
    rgi_df, genome_ids = get_rgi_df(args.rgi_out)

    # filter to whichever criteria you want e.g. perfect only
    if args.perfect_only:
        rgi_df = rgi_df[rgi_df['Cut_Off'] == 'Perfect']

    # encode the TSV

    rgi_encoded = encode_rgi(rgi_df, genome_ids)
    # print(">>>" , ast_df.index)
    try:
        rgi_encoded = rgi_encoded.loc[ast_df.index]
    except Exception as e:
        print("1")
        raise e

    try:
        print("Saving encoded RGI results as rgi_encoded.pkl")
        rgi_encoded[np.isnan(rgi_encoded)] = 0
        rgi_encoded= rgi_encoded.astype(int)
        rgi_encoded.to_pickle('rgi_encoded.pkl')
    except Exception as e:
        print("2")
        raise e
    exit()


    # get alternative encodings
    drug_class_encoded = get_tag_encoded_df(rgi_df, genome_ids, 'Drug Class')
    print("Saving drug class encoded results as drug_class_encoded.pkl")
    drug_class_encoded = drug_class_encoded.loc[ast_df.index]
    drug_class_encoded.to_pickle('drug_class_encoded.pkl')

    amr_family_encoded = get_tag_encoded_df(rgi_df, genome_ids, 'AMR Gene Family')
    amr_family_encoded = amr_family_encoded.loc[ast_df.index]
    print("Saving AMR family encoded results as amr_family_encoded.pkl")
    amr_family_encoded.to_pickle("amr_family_encoded.pkl")



