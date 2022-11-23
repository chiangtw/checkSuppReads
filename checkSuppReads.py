#! /usr/bin/env python

"""
# Usage

checkSuppReads.py NCL_events.tsv file_list.tsv out/

NCL_events.tsv:
    (NCL events) 6-columns TSV: (chr_d, pos_d, strand_d, chr_a, pos_a, strand_a)
    ((circRNA events) 4-columns TSV: (chr, pos1, pos2, strand))

file_list.tsv:
    3-columns TSV: (sample_id, Fastq_1_path, Fastq_2_path)


out
  |-- pseudo_ref/
  |-- sample_1/
  |-- sample_2/
  ...
  |-- sample_N/
  |-- summary.txt

# Dependency

 - python==3
 - bedtools==2.29.0
 - samtools
 - bwa

# Steps

 1. generate pseudo-references
 2. create bwa index
 3. bwa mapping
 4. 


"""

import argparse
import os
import os.path
import io
import subprocess as sp
import csv
import tempfile as tp
import logging
from contextlib import contextmanager
from collections import namedtuple


logging.basicConfig(
    format="{asctime} - {message}",
    level=logging.INFO,
    style='{'
)


class NCLevent:
    JuncSite = namedtuple('JuncSite', ('chr_', 'pos', 'strand'))

    def __init__(self, NCL_data):
        self._data = NCL_data
        self.donor = self.JuncSite(NCL_data[0], int(NCL_data[1]), NCL_data[2])
        self.acceptor = self.JuncSite(NCL_data[3], int(NCL_data[4]), NCL_data[5])

        self.id = "|".join(NCL_data[:6])

    def get_donor_region(self, dist):
        chr_, pos, strand = self.donor

        if strand == '+':
            donor_region = (chr_, pos - dist + 1, pos, strand)
        elif strand == '-':
            donor_region = (chr_, pos, pos + dist - 1, strand)

        return donor_region

    def get_acceptor_region(self, dist):
        chr_, pos, strand = self.acceptor

        if strand == '+':
            acceptor_region = (chr_, pos, pos + dist - 1, strand)
        elif strand == '-':
            acceptor_region = (chr_, pos - dist + 1, pos, strand)

        return acceptor_region


class Bed6:
    def __init__(self, region, name):
        self.region = region
        self.name = name
        self._bed6 = self._to_bed6(region, name)

    @staticmethod
    def _to_bed6(region, name):
        chr_, pos1, pos2, strand = region
        bed6 = (chr_, int(pos1) - 1, int(pos2), name, '.', strand)
        return bed6

    def __str__(self):
        return '\t'.join(map(str, self._bed6))


def get_fasta(genome_file, bed_file, bedtools_bin='bedtools'):
    with tp.NamedTemporaryFile(dir='.') as tmp_file:
        cmd = [bedtools_bin, 'getfasta']
        cmd += ['-fi', genome_file]
        cmd += ['-bed', bed_file]
        cmd += ['-fo', tmp_file.name]
        cmd += ['-name', '-s', '-tab']

        sp.run(cmd)

        with open(tmp_file.name) as fa_in:
            for line in fa_in:
                name, fa_seq = line.rstrip('\n').split('\t')
                yield name, fa_seq


def create_bwa_index(fasta_file, bwa_bin='bwa'):
    cmd = [bwa_bin, 'index', fasta_file]
    result = sp.run(cmd)
    return result


@contextmanager
def cwd(path):
    origin_pwd = os.getcwd()
    os.chdir(path)
    yield
    os.chdir(origin_pwd)


def generate_pseudo_references(NCL_events, genome_file, out_dir, dist=100):
    genome_file = os.path.abspath(genome_file)
    out_dir = os.path.abspath(out_dir)
    pseudo_refs_dir = os.path.join(out_dir, 'pseudo_refs')
    os.makedirs(pseudo_refs_dir, exist_ok=True)

    NCL_file = 'NCL_events.tsv'
    NCL_bed = 'NCL_events.near_junction_region.bed'
    NCL_fasta = 'NCL_events.near_junction_region.fa'

    with cwd(pseudo_refs_dir):
        with open(NCL_file, 'w') as ncl_out, \
                open(NCL_bed, 'w') as bed_out:

            for ncl_ev in NCL_events:
                print(*ncl_ev.donor, *ncl_ev.acceptor, sep='\t', file=ncl_out)

                donor_region = ncl_ev.get_donor_region(dist)
                acceptor_region = ncl_ev.get_acceptor_region(dist)

                donor_bed6 = Bed6(donor_region, f'{ncl_ev.id}_1')
                acceptor_bed6 = Bed6(acceptor_region, f'{ncl_ev.id}_2')

                print(donor_bed6, file=bed_out)
                print(acceptor_bed6, file=bed_out)

        fasta = get_fasta(genome_file, NCL_bed)
        fasta_dict = {name: fa_seq for name, fa_seq in fasta}
        with open(NCL_fasta, 'w') as fa_out:
            for ncl_ev in NCL_events:
                donor_seq = fasta_dict[f'{ncl_ev.id}_1']
                acceptor_seq = fasta_dict[f'{ncl_ev.id}_2']
                
                print(f'>{ncl_ev.id}', file=fa_out)
                print(f'{donor_seq}{acceptor_seq}', file=fa_out)

        create_bwa_index(NCL_fasta)

    return pseudo_refs_dir


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
    parser.add_argument('out_dir')
    parser.add_argument('--index')
    parser.add_argument('-g', '--genome')
    parser.add_argument('-d', '--dist', type=int, default=100)
    parser.add_argument('-t', '--threads', type=int)

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    reader = csv.reader(args.NCL_events, delimiter='\t')
    NCL_events = [NCLevent(data) for data in reader]

    ref_dir = args.index
    if not ref_dir:
        ref_dir = generate_pseudo_references(NCL_events, args.genome, args.out_dir, args.dist)

















