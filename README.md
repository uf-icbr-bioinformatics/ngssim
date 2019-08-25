# ngssim
A suite of simulators for NGS data

## Overview
**ngssim** is a set of four programs that perform different types of
simulations on NGS data. The programs are:

Program               | Opt | Description
----------------------|-----|------------
[simcov](#simcov)     | -C  | Simulate short-read coverage
[simseq](#simseq)     | -S  | Generate a random reference sequence
[simreads](#simreads) | -R  | Generate simulated short reads
[simfastq](#simfastq) | -F  | Generate simulated short reads with quality scores from an existing fastq

The four programs are actually a single executable (simcov.py); the different modes
can be selected using the command-line option in the Opt column, or creating symbolic
links to simcov.py with the appropriate name (only the basename of the link is used).

For example, to invoke simseq you can either use:

```bash
simcov.py -S [options...]
```

or create a symlink having `simseq` as its basename:

```bash
ln -s simcov.py simseq.py
```

and call simseq.py directly.

## Installation

This program is distributed as a single Python2 script, simply save it to a location
in your PATH. Its dependencies are:

* numpy
* The [BIutils](https://github.com/uf-icbr-bioinformatics/BIutils) package.

## simcov

This program places simulated short reads on a genome and computes
the resulting average coverage and other statistics. It can optionally
simulate the presence of variant sites and compute average coverage
over them (useful to plan variant detection experiments).

Options:

Option | Description
-------|------------
  -d   | Run with default paramers (demo).
  -l L | Set read length to L (defalt: 150).
  -n N | Set number of reads to N (default: 10000000).
  -g G | Set target genome size to G (default: 3100000000).
  -p   | Enable paired-end mode (default: single-end).
  -t T | Set insert size to T in paired-end mode (default: 400).
  -f F | Set site frequency to F (default: no site analysis).
  -s S | Set simulation vector size to S (default: 1000000).

Values for -l, -n, -g and -t can be followed by G (for billion) or M (for million).
The value for -f can be expressed as a float or a fraction (e.g. 1/8)

The program simulates N reads of length L, covering a genome of
size G. If -p is specified, reads are simulated as paired-end,
with an insert size of T. Default parameters simulate 10 million reads
on a human-sized genome, and produces output similar to the following:

```
=== Configuration ===
Genome size: 3,100,000,000bp
Number of reads: 10,000,000
Read length: 150bp
Mode: unpaired
Expected average coverage: 0.5X

=== Results ===
Effective average coverage: 0.5X
Max coverage: 5
Bases covered at   1X: 1,182,758,500bp (38.15%)
Bases covered at   5X: 520,800bp (0.02%)
Bases covered at  10X: 0bp (0.00%)
Bases covered at  20X: 0bp (0.00%)
Bases covered at  30X: 0bp (0.00%)
Bases covered at  50X: 0bp (0.00%)
Bases covered at 100X: 0bp (0.00%)
25th percentile - bases covered at 1X: 1,182,758,500bp (38.2%)
50th percentile - bases covered at 2X: 270,878,000bp (8.7%)
75th percentile - bases covered at 4X: 4,898,000bp (0.2%)
```

The simulation is performed on a sequence of one million basepairs
and the results are scaled up to the requested genome size. The size
of the simulated sequence can be changed with the -s option; a larger
size will give more accurate results but will require more time and
memory.

The -f option can be used to simulate the presence of variants (SNPs)
in the reads. Its value is the site frequency, e.g. `-f 1/1000`
would result in one SNP every 1000bp on average. In this
case, the program prints the following additional output:

```
=== SNPs ===
Simulated: 1000
Detected: 378 (37.8%)
  at 1X: 378 (37.8%)
  at 5X: 0 (0.0%)
  at 10X: 0 (0.0%)
  at 20X: 0 (0.0%)
  at 30X: 0 (0.0%)
  at 50X: 0 (0.0%)
  at 100X: 0 (0.0%)
Average squared error: 0.155680669631
```

For each simulated SNP, the program computes the observed allelic frequencies, 
and compares them with the true ones. The last line reports the average squared error in the estimation of the
allelic frequencies of the simulated SNPs. Higher coverage will
result in more accurate frequency estimation and therefore a lower
average.

Note that the number of SNPs is computed on the basis of the size
of the simulation sequence, and not on the true genome size. This has no effect on the displayed percentages and 
average error.

## simseq

This program writes a random sequence to a specified file in FASTA format.

Usage:

```
simseq.py [options] filename.fa
```

where `filename.fa` is the name of the output file, and options are:

Option | Description
-------|------------
 -sn S | Name of sequence (default: seq).
 -l L  | Set sequence length to L (default: 1000000)

The value for -l can be followed by G (for billion) or M (for million).

## simreads

This program extracts short reads from a specified reference sequence
in either single-end or paired-end mode. Quality scores are randomly
generated according to a realistic error model. Output is written to
one or two compressed fastq files.

Usage:

```
simreads.py [options] filename.fa
```

where filename.fa is a FASTA file containing a reference sequence, and
options are:

Option  | Description
--------|------------
  -sn N | Base name for reads (default: seq).
  -nr R | Number of reads (default: 1000000).
  -rl L | Read length (default: 150).
  -i  I | Average insert size (default: 400).
  -is S | Standard deviation of insert size (default: 10).
  -qs S | Average quality at first base (default: 40).
  -qe E | Average quality at last base (default: 30).
  -qv V | Quality standard deviation at last base (default: 10).
  -o  O | Use O as base name for output files (default: reads).
  -p    | Enable paired-end mode.
  -s P  | Simulate presence of P SNPs.
  -so S | Write SNPs to this file (default: None).

The value for -nr can be followed by G (for billion) or M (for million).

Output is written to outfile.fastq.gz, where `outfile` can be set with
the -o option. In paired-end mode (activated with -p) output is written
to outfile_R1.fastq.gz and outfile_R2.fastq.gz.

In paired-end mode, the program simulates the insert size using the average
size specified with -i and a standard deviation specified with -is.

Quality scores are simulated using an exponential decay model where mean
quality scores range from the value set with -qs (at the first base) to the
value set with -qs (at the last base), and the standard deviation ranges from
1 to the value set with -qv (at the last base). For example, using default parameters,
the last base in a read will have a quality score generated from a normal distribution
with mean 30 and standard deviation 10. Bases are subject to random
change according to the simulated quality scores. For example, a base with
a quality score of 20 has a 0.01 probability of being changed.

If -s is used, the program inserts the specified number of SNPs in the reference
sequence before generating the reads. The simulated SNPs can be written to a
file using the -so option. This allows comparing the results of SNP calling
on the simulated reads with the known SNP positions and allelic frequencies, to
evaluate the performance of the SNP caller used.

# simfastq

This program writes randomly-generated reads from a reference sequence, taking
read lengths and qualities from existing fastq file(s). This is useful to generate
sets of random reads for testing programs without altering their quality scores.

Usage:

```
simfastq.py [options] filename.fa infile.fastq[.gz]
```

where `filename.fa` is a reference sequence file in FASTA format, `infile.fastq`
is an (optionally compressed) fastq file, and options are:

Option  | Description
--------|------------
  -o O  | Use O as base name for output file (default: reads.fastq.gz).

For each read in infile, the program will generate a random read from the reference
sequence having the same read name and length, and the same quality scores. Bases
in the read will be mutated at random based on the quality score. For example, if in
the original read a base has a Q score of 20, that base will have a 0.01 probability
of being changed in the output.

If input is paired-end, simply run this program twice on the R1 and R2 files. Output 
is written to file `outfile`.fastq.gz, where outfile can be changed with the -o option.

NOTE: if the reference file contains multiple sequences, there is a small chance
that the output will contain characters from the sequence names. To avoid this,
please index the reference file with samtools:

```
samtools faidx reference.fa
```

If the .fai file for the reference file is present, simfastq will read it automatically 
and will use it to skip the sequence headers.

## Credits
**ngssim** is (c) 2019, A. Riva, ICBR Bioinformatics Core, University of Florida.
