## Usage

```
usage: checkSuppReads.py [-h] [--index INDEX] [-g GENOME] [-d DIST] [-D DIST_FILE] [-l CROSS_LEN] [-m MAP_LEN] [-s SIMILARITY] [-t THREADS] NCL_events file_list out_dir

Checking if there are junction supporting reads for the NCL junction.

positional arguments:
  NCL_events            6-columns TSV: (chr_d, pos_d, strand_d, chr_a, pos_a, strand_a)
  file_list             The file list of samples, consists of 3 columns: (sample_id, path_to_fastq_1, path_to_fastq_2)
  out_dir

optional arguments:
  -h, --help            show this help message and exit
  --index INDEX         Path to the pre-build index, e.g. "./out_dir/pseudo_ref" (default: None)
  -g GENOME, --genome GENOME
  -d DIST, --dist DIST  The extended distance from NCL junction to upstream/downstream. (default: 100)
  -D DIST_FILE, --dist_file DIST_FILE
                        The file of the extended distances. (default: None)
  -l CROSS_LEN, --cross_len CROSS_LEN
                        The minimal length of bases across the NCL junction. (default: 10)
  -m MAP_LEN, --map_len MAP_LEN
                        . (default: 20)
  -s SIMILARITY, --similarity SIMILARITY
                        . (default: 0.8)
  -t THREADS, --threads THREADS
                        . (default: 1)

```

#### Format of DIST_FILE

The DIST_FILE is consist of three columns: (NCL_id, dist_donor, dist_acceptor).

eg.
```
chr9|137145902|+|chr9|136362180|-    100  100
chr5|149845914|+|chr5|149832656|+    100  100
chr10|133381659|+|chr10|133379999|+  100  100
chrX|107088436|-|chrX|107115556|-    100  100
chr1|28476291|+|chr1|28466382|+      100  100
```

See [examples](examples).

## Requirements

- python==3
- bedtools
- samtools
- bwa
