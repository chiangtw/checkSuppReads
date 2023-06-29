#! /usr/bin/env python


import argparse
import csv
import os
import os.path
from contextlib import contextmanager
from itertools import groupby
from operator import itemgetter


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


class PseudoRefDistDB:
    def __init__(self, dist=100, dist_file=None):
        self._dist = dist
        self._dist_file = dist_file

        self._dist_db = {}

        if self._dist_file:
            self._parse_dist_file(self._dist_file)

    def _parse_dist_file(self, dist_file):
        with open(dist_file) as f_in:
            for line in f_in:
                id_, dist_d, dist_a = line.rstrip('\n').split('\t')
                dist_d = int(dist_d)
                dist_a = int(dist_a)

                self._dist_db[id_] = (dist_d, dist_a)

    def get(self, id_):
        return self._dist_db.get(id_, (self._dist, self._dist))


def filter_junction_reads(junc_reads_file, out_file, dist_db, cross_junc_threshold, map_len_threshold, similarity_threshold):
    with open(junc_reads_file) as f_in, open(out_file, 'w') as out:
        for line in f_in:
            data = line.rstrip('\n').split('\t')

            ref_name = data[1]
            pos = int(data[2])
            Z3 = int(data[3])
            map_len = int(data[4])
            ZS = float(data[7])

            ref_len_d, _ = dist_db.get(ref_name)

            if (map_len >= map_len_threshold) and \
                    (ZS >= similarity_threshold) and \
                    (pos <= ref_len_d - cross_junc_threshold + 1) and \
                    (Z3 >= ref_len_d + cross_junc_threshold):

                print(*data, sep='\t', file=out)

    return os.path.abspath(out_file)


def retain_uniq_read_ref(s1_s2_file, out_file):
    with open(s1_s2_file) as f_in:
        all_data = [line.rstrip('\n').split('\t') for line in f_in]

    all_read_ref_pairs = sorted(set((data[0], data[1]) for data in all_data))

    with open(out_file, 'w') as out:
        for read_id, read_gp in groupby(all_read_ref_pairs, key=itemgetter(0)):
            read_gp = list(read_gp)

            if len(read_gp) == 1:
                print(*read_gp[0], sep='\t', file=out)

    return os.path.abspath(out_file)


def merge_all_supporting_reads(s1_s2_uniq_file, out_file):
    with open(s1_s2_uniq_file) as f_in:
        all_data = [line.rstrip('\n').split('\t') for line in f_in]

    all_ref_read_pairs = sorted([(data[1], data[0]) for data in all_data])
    
    with open(out_file, 'w') as out:
        for ref_id, ref_gp in groupby(all_ref_read_pairs, key=itemgetter(0)):
            ref_gp = list(ref_gp)
            all_reads = [pair[1] for pair in ref_gp]
            print(ref_id, len(all_reads), ','.join(all_reads), sep='\t', file=out)

    return os.path.abspath(out_file)


def modify_conditions(results_dir, sample_id, dist_db, cross_len, map_len, similarity):
    sample_dir = os.path.join(results_dir, sample_id)

    with cwd(sample_dir):
        s1_s2_filtered_file = filter_junction_reads('all_junc_reads_s1_s2.tsv', 'all_junc_reads_s1_s2.filtered.tsv', dist_db, cross_len, map_len, similarity)
        s1_s2_uniq_file = retain_uniq_read_ref(s1_s2_filtered_file, 'all_junc_reads_s1_s2.tsv.uniq_read_ref')
        read_count_file = merge_all_supporting_reads(s1_s2_uniq_file, 'all_junc_reads.count')


def create_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'file_list',
        type=argparse.FileType('r'),
        help=('The file list of samples, consists of 3 columns: '
              '(sample_id, path_to_fastq_1, path_to_fastq_2)')
    )
    parser.add_argument('results_dir', help='The results directory from "checkSuppReads.py".')
    parser.add_argument('index_dir', help='Path to the pre-build index, e.g. "./out_dir/pseudo_ref"')
    parser.add_argument('-l', '--cross_len', type=int, default=10, help='The minimal length of bases across the NCL junction.')
    parser.add_argument('-m', '--map_len', type=int, default=20, help='.')
    parser.add_argument('-s', '--similarity', type=float, default=0.8, help='.')

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    dist_db = PseudoRefDistDB(dist_file=os.path.join(args.index_dir, 'NCL_events.dist.tsv'))

    sample_ids = get_sample_ids(args.file_list)

    for sample_id in sample_ids:
        modify_conditions(args.results_dir, sample_id, dist_db, args.cross_len, args.map_len, args.similarity)
