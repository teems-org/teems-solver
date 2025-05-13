import numpy as np
import pandas as pd
from glob import glob
from toolz import compose
from os.path import join
from os import makedirs

read_int_array = compose(np.array, lambda x: list(map(int, x)))


def extract_name(filename):
    with open(filename) as sample_file:
        name = next(sample_file).strip().split(' ')[-2]
    return name


def split_variables(strip_list):

    if len(strip_list) == 0:
        return []

    cur_header = strip_list[0]
    n_to_expect = int(cur_header.split(' ')[0])
    cur_name = cur_header.split(' ')[1].split(':')[0]

    # Go past header
    strip_list = strip_list[1:]

    cur_vars = strip_list[:n_to_expect]

    return [(cur_name, cur_vars)] + split_variables(
        strip_list[n_to_expect:])


def parse_vector(vec_str, typename=int):
    return np.array([typename(x) for x in vec_str.split(' ')])


def parse_array(vec_str_list, typename=int):
    return np.stack([parse_vector(x) for x in vec_str_list])


def read_mds_file(mds_filename):

    # The [1:] is to skip the header
    entries = [x.strip() for x in open(mds_filename)][1:]

    nsetspace, nvar, nvarele, nset = map(int, entries)

    return nsetspace, nvar, nvarele, nset


def read_bin_file(bin_filename):

    bin_file = [x.strip() for x in open(bin_filename)][1:]
    xc = np.array([float(x) for x in bin_file])

    return xc


def read_var_file(var_filename):

    var_file = [x.strip() for x in open(var_filename)][1:]

    split_var_file = split_variables(var_file)
    split_var_file = {x: y for (x, y) in split_var_file}

    # Further split the arrays:
    split_var_file['setid'] = parse_array(split_var_file['setid'])
    split_var_file['antidims'] = parse_array(split_var_file['antidims'])

    # And parse the rest into the right formats
    int_arrays = ['begadd', 'size', 'matsize', 'level_par', 'change_real',
                  'suplval', 'gltype']

    for cur_to_convert in int_arrays:
        split_var_file[cur_to_convert] = np.array(
            [int(x) for x in split_var_file[cur_to_convert]])

    split_var_file['glval'] = np.array([float(x) for x in
                                        split_var_file['glval']])
    split_var_file['cofname'] = np.array(split_var_file['cofname'])

    return split_var_file


def read_sel_file(sel_filename):

    sel_file = [x.strip() for x in open(sel_filename)][1:]

    split_sel_file = split_variables(sel_file)
    split_sel_file = {x: y for (x, y) in split_sel_file}

    split_sel_file['setele'] = np.array(split_sel_file['setele'])
    split_sel_file['setsh'] = parse_array(split_sel_file['setsh'])

    return split_sel_file


def read_set_file(set_filename):

    set_file = [x.strip() for x in open(set_filename)][1:]

    split_set_file = split_variables(set_file)
    split_set_file = {x: y for (x, y) in split_set_file}

    funs_to_apply = {
        'header': np.array,
        'fileid': read_int_array,
        'setname': np.array,
        'readele': np.array,
        'begadd': read_int_array,
        'size': read_int_array,
        'subsetid': parse_array,
        'intertemp': read_int_array,
        'intsup': read_int_array,
        'regional': read_int_array,
        'regsup': read_int_array
    }

    split_set_file = {x: funs_to_apply[x](y) for x, y in
                      split_set_file.items()}

    return split_set_file


def save_individual_csvs(target_dir, make_dir_if_required=True, **params):

    if make_dir_if_required:
        makedirs(target_dir, exist_ok=True)

    for cur_name, cur_value in params.items():

        pd.DataFrame(cur_value).to_csv(join(target_dir, cur_name + '.csv'))


if __name__ == '__main__':

    import sys

    base_dir = sys.argv[1]
    target_dir = sys.argv[2]

    makedirs(target_dir, exist_ok=True)

    all_files = glob(join(base_dir, '*'))

    lookup = {extract_name(x): x for x in all_files}

    nsetspace, nvar, nvarele, nset = read_mds_file(lookup['.mds'])
    xc = read_bin_file(lookup['.bin'])
    split_var_file = read_var_file(lookup['.var'])
    split_sel_file = read_sel_file(lookup['.sel'])
    split_set_file = read_set_file(lookup['.set'])

    # Also save csv versions
    for cur_name, cur_dict in zip(['var_csvs', 'sel_csvs', 'set_csvs'],
                                  [split_var_file, split_sel_file,
                                   split_set_file]):

        save_individual_csvs(join(target_dir, f'{cur_name}'),
                             **cur_dict)

    # Save the other ones individually
    save_individual_csvs(join(target_dir, 'bin_csvs'), xc=xc)
    save_individual_csvs(join(target_dir, 'mds_csvs'), nsetspace=[nsetspace],
                         nvar=[nvar], nvarele=[nvarele], nset=[nset])
