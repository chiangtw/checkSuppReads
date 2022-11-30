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
 4. filtering


"""

import argparse
import os
import os.path
import io
import subprocess as sp
import csv
import tempfile as tp
import re
import logging
from contextlib import contextmanager
from collections import namedtuple, defaultdict
from operator import itemgetter
from itertools import chain, groupby


logging.basicConfig(
    format="{asctime} - {message}",
    level=logging.INFO,
    style='{'
)


BEDTOOLS_BIN = 'bedtools'
SAMTOOLS_BIN = 'samtools'
BWA_BIN = 'bwa'


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


@contextmanager
def cwd(path):
    origin_pwd = os.getcwd()
    os.chdir(path)
    yield
    os.chdir(origin_pwd)


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


def get_fasta(genome_file, bed_file, bedtools_bin=BEDTOOLS_BIN):
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


def create_bwa_index(fasta_file, bwa_bin=BWA_BIN):
    cmd = [bwa_bin, 'index', fasta_file]
    result = sp.run(cmd)
    return result


def generate_pseudo_references(NCL_events, genome_file, out_dir, dist=100):
    logging.info('Generating pseudo-references')
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

    return os.path.join(pseudo_refs_dir, NCL_fasta)


def bwa_mapping(index_file, fastq_file, out_file, threads=1, bwa_bin=BWA_BIN, samtools_bin=SAMTOOLS_BIN):
    logging.info(f'bwa mapping for {fastq_file}')

    cmd_1 = [
        bwa_bin, 'mem',
        '-t', str(threads),
        '-O', '100',
        index_file,
        fastq_file
    ]

    cmd_2 = [samtools_bin, 'view', '-bh', '-']

    with open(out_file, 'wb') as bam_out:
        p1 = sp.Popen(cmd_1, stdout=sp.PIPE, encoding='utf-8')
        p2 = sp.Popen(cmd_2, stdin=p1.stdout, stdout=bam_out)

        p1.wait()
        p2.wait()

    return os.path.abspath(out_file)


class CIGAR:
    """Use to parse the CIGAR string."""

    _CIGAR_PAT = re.compile(r'(([0-9]+)([MIDNSHP=X]))')

    def __init__(self, cigar_str):
        self._cigar_str = cigar_str
        self._parse()

    def _parse(self):

        if self._cigar_str == '*':
            self._cigar_list = []
            self._cigar_dict = {}
        else:
            parsed_result = re.findall(self._CIGAR_PAT, self._cigar_str)

            self._cigar_list = [cigar for cigar, _, _ in parsed_result]

            self._cigar_dict = defaultdict(list)
            for _, num, op in parsed_result:
                self._cigar_dict[op].append(int(num))

    def __getitem__(self, key):
        if isinstance(key, int):
            return self._cigar_list[key]

        elif isinstance(key, str):
            return self._cigar_dict.get(key, [])

        else:
            pass

    def __repr__(self):
        return str(self._cigar_list)

    def __str__(self):
        return self._cigar_str


class SamFormat:
    _optional_field_pattern = re.compile(r'([A-Za-z][A-Za-z0-9]):([AifZHB]):(.+)')
    _optional_field_tuple = namedtuple('OptField', ['tag', 'type_', 'value'])

    def __init__(self, sam_string):
        self._sam_string = sam_string.rstrip('\n')
        self._init()
        self._parse()

    def _init(self):
        self.is_header = False

        self.qname = None
        self.flag = None
        self.rname = None
        self.pos = None
        self.mapq = None
        self.cigar = None
        self.rnext = None
        self.pnext = None
        self.tlen = None
        self.seq = None
        self.qual = None
        self.optional_fields = None

    def _parse(self):
        if self._sam_string.startswith('@'):
            self.is_header = True
        else:
            data = self._sam_string.split('\t')

            self.qname = data[0]
            self.flag = int(data[1])
            self.rname = data[2]
            self.pos = int(data[3])
            self.mapq = int(data[4])
            self.cigar = CIGAR(data[5])
            self.rnext = data[6]
            self.pnext = int(data[7])
            self.tlen = int(data[8])
            self.seq = data[9]
            self.qual = data[10]

            self.optional_fields = self._parse_optional_fields(data[11:])

    def _parse_optional_fields(self, fields):
        fields_dict = {}
        for field in fields:
            m = re.search(self._optional_field_pattern, field)
            if m:
                tag, type_, value = m.groups()
                fields_dict[tag] = self._optional_field_tuple(tag, type_, value)

        return fields_dict

    def __repr__(self):
        return self._sam_string

    def __str__(self):
        return self._sam_string

    @property
    def is_unmapped(self):
        return str(self.cigar) == '*'


def generate_Z3_tag(sam_data):
    ref_consumes_items = itemgetter('M', 'D', 'N', '=', 'X')(sam_data.cigar)
    ref_consumes = sum(chain.from_iterable(ref_consumes_items))
    Z3 = sam_data.pos + ref_consumes - 1
    return "Z3:i:{}".format(Z3)


def calc_similarity(cigar, MD):
    total_matches = sum(cigar['M'])
    perfect_matches = sum(map(int, re.findall(r'([0-9]+)', MD)))

    similarity = perfect_matches / total_matches

    return similarity


def generate_ZS_tag(sam_data):
    similarity = calc_similarity(
        sam_data.cigar,
        sam_data.optional_fields['MD'].value
    )
    return "ZS:f:{}".format(round(similarity, 4))


@contextmanager
def bam_reader(bam_file, samtools_bin=SAMTOOLS_BIN):
    cmd = [samtools_bin, 'view', '-h', bam_file]

    with sp.Popen(cmd, stdout=sp.PIPE, encoding='utf-8') as p:
        yield iter(p.stdout.readline, '')


@contextmanager
def bam_writer(out_file, samtools_bin=SAMTOOLS_BIN):
    cmd = [samtools_bin, 'view', '-bh', '-o', out_file, '-']

    with sp.Popen(cmd, stdin=sp.PIPE, encoding='utf-8') as p:
        yield p.stdin


def append_Z3_ZS_tag(bam_file, out_file):
    logging.info(f'appending Z3 & ZS tag to {bam_file}')
    with bam_reader(bam_file) as reader:
        with bam_writer(out_file) as writer:
            for line in reader:
                sam_data = SamFormat(line)

                if sam_data.is_header:
                    print(sam_data, file=writer)
                else:
                    if sam_data.is_unmapped:
                        print(sam_data, file=writer)
                    else:
                        Z3_tag = generate_Z3_tag(sam_data)
                        ZS_tag = generate_ZS_tag(sam_data)
                        print(sam_data, Z3_tag, ZS_tag, sep='\t', file=writer)

    return os.path.abspath(out_file)


def get_uniq_matches(bam_file, out_file):
    logging.info('getting unique matches')
    with bam_reader(bam_file) as reader:
        with bam_writer(out_file) as writer:

            qname = None
            sam_gp = []

            for line in reader:
                sam_data = SamFormat(line)

                if sam_data.is_header:
                    print(sam_data, file=writer)
                else:
                    if qname is None:
                        qname = sam_data.qname
                        sam_gp.append(sam_data)
                    else:
                        if sam_data.qname == qname:
                            sam_gp.append(sam_data)
                        else:
                            flags = set(sam.flag for sam in sam_gp)
                            for flag in flags:
                                if flag & 2048:
                                    break
                            else:
                                for sam in sam_gp:
                                    print(sam, file=writer)

                            # init for next group
                            qname = sam_data.qname
                            sam_gp = [sam_data]

    return os.path.abspath(out_file)


def get_junc_reads(bam_file, out_file, ref_len, cross_junc_threshold, map_len_threshold, similarity_threshold):
    logging.info('getting junction reads')
    with bam_reader(bam_file) as reader, open(out_file, 'w') as out:
        for line in reader:
            sam_data = SamFormat(line)

            if not sam_data.is_header:
                if not sam_data.is_unmapped:
                    pos = sam_data.pos
                    Z3 = int(sam_data.optional_fields['Z3'].value)
                    ZS = float(sam_data.optional_fields['ZS'].value)

                    map_len = Z3 - sam_data.pos + 1

                    if (map_len >= map_len_threshold) and \
                            (ZS >= similarity_threshold) and \
                            (pos <= ref_len - cross_junc_threshold + 1) and \
                            (Z3 >= ref_len + cross_junc_threshold):

                        print(
                            sam_data.qname,
                            sam_data.rname,
                            sam_data.pos,
                            Z3,
                            map_len,
                            sam_data.cigar,
                            sam_data.optional_fields['MD'].value,
                            ZS,
                            sep='\t',
                            file=out
                        )

    return os.path.abspath(out_file)


def merge_junction_reads(fastq1_junc_reads, fastq2_junc_reads, out_file):
    all_junc_reads = []

    with open(fastq1_junc_reads) as f_in:
        for line in f_in:
            data = line.rstrip('\n').split('\t')
            all_junc_reads.append(data + ['1'])

    with open(fastq2_junc_reads) as f_in:
        for line in f_in:
            data = line.rstrip('\n').split('\t')
            all_junc_reads.append(data + ['2'])

    all_junc_reads = sorted(all_junc_reads, key=itemgetter(0, 1, 8))

    with open(out_file, 'w') as out:
        for data in all_junc_reads:
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


def check_supporting_reads(index_file, sample_id, fastq1, fastq2, out_dir, threads=1, dist=100):
    logging.info('Checking supporting reads')
    index_file = os.path.abspath(index_file)
    fastq1 = os.path.abspath(fastq1)
    fastq2 = os.path.abspath(fastq2)

    sample_dir = os.path.join(out_dir, sample_id)
    os.makedirs(sample_dir, exist_ok=True)

    fastq1_dir = os.path.join(sample_dir, 'fastq1')
    os.makedirs(fastq1_dir, exist_ok=True)

    fastq2_dir = os.path.join(sample_dir, 'fastq2')
    os.makedirs(fastq2_dir, exist_ok=True)

    with cwd(fastq1_dir):
        logging.info(f'Checking fastq1 of {sample_id}')
        fastq1_bam = bwa_mapping(index_file, fastq1, 'Aligned.out.bam', threads)
        fastq1_Z3_ZS_bam = append_Z3_ZS_tag(fastq1_bam, 'Aligned.out.Z3.ZS.bam')
        fastq1_uniq_bam = get_uniq_matches(fastq1_Z3_ZS_bam, 'Aligned.out.Z3.ZS.uniq_matches.bam')
        fastq1_junc_reads = get_junc_reads(fastq1_uniq_bam, 'Aligned.out.Z3.ZS.uniq_matches.bam.junc_reads', dist, 10, 20, 0.8)

    with cwd(fastq2_dir):
        logging.info(f'Checking fastq2 of {sample_id}')
        fastq2_bam = bwa_mapping(index_file, fastq2, 'Aligned.out.bam', threads)
        fastq2_Z3_ZS_bam = append_Z3_ZS_tag(fastq2_bam, 'Aligned.out.Z3.ZS.bam')
        fastq2_uniq_bam = get_uniq_matches(fastq2_Z3_ZS_bam, 'Aligned.out.Z3.ZS.uniq_matches.bam')
        fastq2_junc_reads = get_junc_reads(fastq2_uniq_bam, 'Aligned.out.Z3.ZS.uniq_matches.bam.junc_reads', dist, 10, 20, 0.8)

    with cwd(sample_dir):
        logging.info(f'Merging results for {sample_id}')
        s1_s2_file = merge_junction_reads(fastq1_junc_reads, fastq2_junc_reads, 'all_junc_reads_s1_s2.tsv')
        s1_s2_uniq_file = retain_uniq_read_ref(s1_s2_file, 'all_junc_reads_s1_s2.tsv.uniq_read_ref')
        read_count_file = merge_all_supporting_reads(s1_s2_uniq_file, 'all_junc_reads.count')

    return read_count_file


def output_summary(all_results, out_dir):
    summary = defaultdict(dict)
    with cwd(out_dir):
        for sample_id, result_file in all_results:
            with open(result_file) as f_in:
                for line in f_in:
                    data = line.rstrip('\n').split('\t')
                    NCL_id = data[0]
                    read_count = data[1]

                    summary[NCL_id][sample_id] = read_count

        with open('summary.txt', 'w') as out:
            all_sample_ids = [sample_id for sample_id, _ in all_results]
            print('NCL_id', *all_sample_ids, sep='\t', file=out)

            for NCL_id in summary:
                all_read_count = []
                for sample_id in all_sample_ids:
                    read_count = summary[NCL_id].get(sample_id, '0')
                    all_read_count.append(read_count)

                print(NCL_id, *all_read_count, sep='\t', file=out)


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
    parser.add_argument('-t', '--threads', type=int, default=1)

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    ncl_reader = csv.reader(args.NCL_events, delimiter='\t')
    NCL_events = [NCLevent(data) for data in ncl_reader]

    index_file = args.index
    if not index_file:
        index_file = generate_pseudo_references(NCL_events, args.genome, args.out_dir, args.dist)

    all_results = []
    file_reader = csv.reader(args.file_list, delimiter='\t')
    for sample_id, fastq1, fastq2 in file_reader:
        result_file = check_supporting_reads(index_file, sample_id, fastq1, fastq2, args.out_dir, args.threads, args.dist)
        all_results.append([sample_id, result_file])

    output_summary(all_results, args.out_dir)
