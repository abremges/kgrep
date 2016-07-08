kgrep
=====

Roughly speaking, a 'grep' for k-mers.

The primary use case is the recruitment of metagenomic proxy reads, starting from a set of single cell reads. Assembling these instead of the single cell reads into a so-called metagenomic proxy assembly usually results in a nicer assembly.

You can also use the tool to quickly identify known contaminants or, because it is insanely fast, screen multiple metagenomes for the abundance of any genome of interest.

## Installation
    git clone https://github.com/metagenomics/kgrep.git
    cd kgrep && make

## Usage
    kgrep [options] <seed> <in1.fq> [in2.fq]

       <seed>     recruitment seed (e.g. single-cell reads or reference sequence)
       <in1.fq>   read pool to recruit from (e.g. metagenomic reads)
       [in2.fq]   optionally: #2 mates, if <in1.fq> contains #1 mates

    Options:

       -p         <in1.fq> consists of interleaved mates (i.e. shuffled reads)
       -x         enable on-the-fly recruitment seed expansion (GREEDY and EXPERIMENTAL)

       -k INT     k-mer size (15 <= k <= 21), requires 4^(k-1) bytes of RAM (15:256M, 16:1G, 17:4G, 18:16G, ...) [17]

       -e FLOAT   default: automatically determine the expected error rate based on Phred quality scores with an upper limit [0.01]
       -d INT     alternatively: min. allowed edit distance (q-gram lemma applied: hits_per_read = readlength - k+1 - d*k) [0]
       -n INT     alternatively: fixed number of required hits per read [0]

       -i FILE    ignore sequence patterns in FILE (e.g. adapter sequences) [null]
       -v         select non-matching sequences

======
*Heavily relies on on [klib](https://github.com/attractivechaos/klib).
Source code in part adopted from – or inspired by – [Heng Li](https://github.com/lh3).*
