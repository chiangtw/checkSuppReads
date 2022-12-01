## Usage

```
usage: checkSuppReads.py [-h] [--index INDEX] [-g GENOME] [-d DIST] [-t THREADS] NCL_events file_list out_dir

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
  -t THREADS, --threads THREADS
                        (default: 1)

```

