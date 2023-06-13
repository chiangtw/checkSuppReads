#! /usr/bin/env python


import argparse
import csv
import os
import os.path
from collections import defaultdict
from contextlib import contextmanager


def get_NCL_ids(NCL_events_file):
    ncl_reader = csv.reader(NCL_events_file, delimiter='\t')
    NCL_ids = ['|'.join(data) for data in ncl_reader]
    return NCL_ids


def get_sample_ids(file_list):
    file_reader = csv.reader(args.file_list, delimiter='\t')
    sample_ids = [sample_id for sample_id, _, _ in file_reader]
    return sample_ids


@contextmanager
def cwd(path):
    origin_pwd = os.getcwd()
    os.chdir(path)
    yield
    os.chdir(origin_pwd)


def get_all_results(sample_ids, results_dir):
    all_results = defaultdict(dict)
    for sample_id in sample_ids:
        result_file = os.path.join(results_dir, sample_id, 'all_junc_reads.count')
        with open(result_file) as f_in:
            for line in f_in:
                data = line.rstrip('\n').split('\t')
                NCL_id = data[0]
                read_count = data[1]

                all_results[NCL_id][sample_id] = read_count

    return all_results


def output_summary(NCL_ids, sample_ids, all_results, out_file):
    print('NCL_id', *sample_ids, sep='\t', file=out_file)

    for ncl_id in NCL_ids:
        all_read_count = []
        for sample_id in sample_ids:
            read_count = all_results[ncl_id].get(sample_id, '0')
            all_read_count.append(read_count)

        print(ncl_id, *all_read_count, sep='\t', file=out_file)


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'NCL_events',
        type=argparse.FileType('r'),
        help='6-columns TSV: (chr_d, pos_d, strand_d, chr_a, pos_a, strand_a)'
    )
    parser.add_argument(
        'file_list',
        type=argparse.FileType('r'),
        help=('The file list of samples, consists of 3 columns: '
              '(sample_id, path_to_fastq_1, path_to_fastq_2)')
    )
    parser.add_argument('results_dir', help='The results directory from "checkSuppReads.py".')
    parser.add_argument('out_summary_file', type=argparse.FileType('w'))

    return parser


if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()

    sample_ids = get_sample_ids(args.file_list)
    all_results = get_all_results(sample_ids, args.results_dir)

    NCL_ids = get_NCL_ids(args.NCL_events)

    output_summary(NCL_ids, sample_ids, all_results, args.out_summary_file)
