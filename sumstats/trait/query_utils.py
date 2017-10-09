import argparse
import numpy as np
import sumstats.utils.utils as utils


TO_LOAD_DSET_HEADERS = ['snp', 'pval', 'chr', 'or', 'bp', 'effect', 'other']
TO_STORE_DSETS = ['snp', 'pval', 'chr', 'or', 'bp', 'effect', 'other']
SNP_DSET = 'snp'
BP_DSET = 'bp'
PVAL_DSET = 'pval'

def get_dsets_from_file(f, dsets):

    # initialize dictionary of datasets
    dict_of_dsets = {dset : [] for dset in dsets}

    for trait, trait_group in f.items():
        dict_trait_dsets = get_dsets_from_trait_group(trait_group, dsets)
        for dset in dict_of_dsets:
            dict_of_dsets[dset].extend(dict_trait_dsets[dset])

    for dset in dsets:
        dict_of_dsets[dset] = np.array(dict_of_dsets[dset])

    return dict_of_dsets


def get_dsets_from_trait_group(trait_group, dsets):

    # initialize empty array for each dataset
    dict_of_dsets = {dset_name : [] for dset_name in dsets}

    for study, study_group in trait_group.items():
        for dset_name in dsets:
            dataset = dict_of_dsets[dset_name]
            dataset.extend(get_dset_from_group(dset_name, study_group, study))

    for dset in dsets:
        dict_of_dsets[dset] = np.array(dict_of_dsets[dset])

    return dict_of_dsets


def get_dset_from_group(dset_name, group, extra_array_filler=None):
    array = utils.get_dset(group, dset_name)
    if (array is None) and (extra_array_filler is not None):
        # pval is never empty
        pval = utils.get_dset(group, PVAL_DSET)
        array = [extra_array_filler for _ in range(len(pval))]
    return array


def argument_checker():
    args = argument_parser()

    if args.query == "1":
        if args.trait is None:
            print("Query 1 -- Retrieve all information for a specific trait -- ")
            print("input: query number (1) and trait name")
            raise SystemExit(1)
    elif args.query == "2":
        if args.study is None or args.trait is None:
            print("Query 2 -- Retrieve all the information for a specific study --")
            print("input: query number (2), study name and trait name")
            raise SystemExit(1)
    else:
        print("Wrong input")
        raise SystemExit(1)


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-h5file', help='The name of the HDF5 file', required=True)
    parser.add_argument('-query', help='The nmber of the query to perform', required=True)
    parser.add_argument('-study', help='The study I am looking for')
    parser.add_argument('-trait', help='The trait I am looking for')
    parser.add_argument('-snp', help='Filter by SNP')
    parser.add_argument('-chr', help='Filter by chromosome')
    parser.add_argument('-pu', help='The upper limit for the p-value')
    parser.add_argument('-pl', help='The lower limit for the p-value')

    return parser.parse_args()