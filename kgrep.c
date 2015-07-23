/* The MIT License

   Copyright (c) 2014, 2015 Andreas Bremges <andreas@cebitec.uni-bielefeld.de>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#define PROGRAM "kgrep"
#define DESCR   "grep with k-mers"
#define VERSION "0.7.0-alpha"
#define AUTHOR  "Andreas Bremges"
#define CONTACT "andreas@cebitec.uni-bielefeld.de" // TODO Too much information...

#include <math.h>     // pow
#include <stdint.h>   // uint64_t
#include <stdio.h>    // printf
#include <sys/time.h> // gettimeofday
#include <unistd.h>   // getopt
#include <zlib.h>     // gzip

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef uint64_t task_t;
#define IGNORE  (task_t) 2
#define WHITE   (task_t) 1
#define RECRUIT (task_t) 0

#define FLAG_P 0x01
#define FLAG_X 0x02
#define FLAG_R 0x04
#define FLAG_D 0x08
#define FLAG_N 0x10
#define FLAG_I 0x20
#define FLAG_C 0x40
#define FLAG_V 0x80
// TODO Replace flag stuff with nice struct

#define LOGMOD 100000

static double start;
static uint64_t *kmers;
static int k = 17;
static int d = 1;
static double max_error = 0.01;
static int flag = 0;
static int phred_offset = 64; // TODO

// A = 00, C = 01, G = 10, T = 11, i.e. XOR with 3 -> compl.
unsigned char seq_fwd_table[128] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

// TODO Taken from ...
double realtime() {
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp, &tzp);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}

// TODO Taken from ...
void stk_printseq(const kseq_t *s) {
    fputc(s->qual.l? '@' : '>', stdout);
    fputs(s->name.s, stdout);
    if (s->comment.l) {
        fputc(' ', stdout);
        fputs(s->comment.s, stdout);
    }
    fputc('\n', stdout);
    fputs(s->seq.s, stdout);
    fputc('\n', stdout);
    if (s->qual.l) {
        fputs("+\n", stdout);
        fputs(s->qual.s, stdout);
        fputc('\n', stdout);
    }
}

int get_min_hits(const kseq_t *seq) {
	if (flag&FLAG_D) return seq->seq.l - k+1 - d*k; // calculate based on min allowed edit distance
	if (flag&FLAG_N) return d; // fixed number of required hits (d overloaded)
	if (seq->qual.l) {
		double p_sum = 0; // TODO init with 1 to round up? Harder to read.
	    for (int i = 0; i < seq->qual.l; ++i) {
	        // Q = Phred quality score, P = base-calling error probability
	        // Then: Q = -10*log_10(P) <=> P = 10^(-Q/10)
	        int q_e = seq->qual.s[i]-phred_offset;
	        double p_e = pow(10,(-q_e/10.0));
	        p_sum += p_e;
	    }
	    //fprintf(stderr, "[%.2f] %f | %f\n", (realtime()-start), seq->seq.l*max_error, p_sum);
	    int exp_max_e = seq->seq.l*max_error+1;
	    int exp_phred_e = p_sum+1; // TODO We calculate the ceiling here rather than taking the double as inout for formula - think about it!
	    // Pick min(auto,max) as the guesstimate
	    if (exp_max_e <= exp_phred_e) {
	    	return seq->seq.l - k+1 - exp_max_e*k;
	    } else {
	    	return seq->seq.l - k+1 - exp_phred_e*k;
	    }
	} else {
		int exp_max_e = seq->seq.l*max_error+1;
		return seq->seq.l - k+1 - exp_max_e*k;
	}
}

int hash_read(const kseq_t *seq, const task_t task) {
	int index = 0;
	int khits = 0;

	uint64_t forward = 0;
	uint64_t reverse = 0;

	for (int i = 0; i < seq->seq.l; ++i) {
		index++;
		int c = seq_fwd_table[(int) seq->seq.s[i]];
		if (c < 4) {
			uint64_t mask = (1ULL<<(k<<1))-1;
			forward = (((forward << 2) | c) & mask); // see below
			//if (mask > 0) forward &= mask; <- only a problem if k = 32
			uint64_t shft = (k<<1)-2;
			reverse = ((reverse >> 2) | (((uint64_t) (c^3)) << shft));
		} else {
			index = 0; //return 0;
		}
		if (index >= k) {
			if (task > 0) { // fill hash
				uint64_t z = ((forward < reverse) ? forward : reverse);
				if ((kmers[z>>5]>>((z&31)<<1)&3) < task) {
					kmers[z>>5] |= (task << ((z&31)<<1));
                    khits++;
				}
			} else { // recruitment step
				uint64_t z = ((forward < reverse) ? forward : reverse);
				int status = (kmers[z>>5]>>((z&31)<<1)&3);
				if (status == WHITE) {
					khits++;
				}
			}
		}
	}
	return khits;
}

// TODO Stop sliding if recruitment criterion is met?
void recruit_pe(const char *file1, const char *file2, const task_t task) {

	gzFile fp1 = strcmp(file1, "-") ? gzopen(file1, "r") : gzdopen(fileno(stdin), "r");
	kseq_t *fqrec1 = kseq_init(fp1);

	gzFile fp2 = NULL;
	kseq_t *fqrec2 = NULL;
	kstream_t *tmp = NULL;

	if (file2) {
        fprintf(stderr, "[%.2f] Now recruiting paired reads from %s and %s\n", (realtime()-start), file1, file2);
		fp2 = gzopen(file2, "r");
		fqrec2 = kseq_init(fp2);
	} else {
        fprintf(stderr, "[%.2f] Now recruiting interleaved pairs from %s\n", (realtime()-start), file1);
		fqrec2 = kseq_init(fp1);
		tmp = fqrec2->f;		// Maybe this can be handled more elegantly?
		fqrec2->f = fqrec1->f;	// fqrec1 and fqrec2 now share one kstream.
	}

	static int counter, counter2;
	while (kseq_read(fqrec1) >= 0 && kseq_read(fqrec2) >= 0) {
		counter += 2;
		int n1, n2;
		if ((n1 = get_min_hits(fqrec1)) > 0 && (n2 = get_min_hits(fqrec2)) > 0) {
			int khits1, khits2;
			if ((khits1 = hash_read(fqrec1, task)) >= n1 && (khits2 = hash_read(fqrec2, task)) >= n2) { // TODO allow OR and AND when recruiting pairs
				if (!(flag&FLAG_V)) {
					counter2 += 2;
					stk_printseq(fqrec1);
					stk_printseq(fqrec2);
				}
				if (flag&FLAG_X) {
					hash_read(fqrec1, WHITE);
					hash_read(fqrec2, WHITE);
				}
			} else if (flag&FLAG_V) {
				counter2 += 2;
				stk_printseq(fqrec1);
				stk_printseq(fqrec2);
			}
		}
		if (counter%LOGMOD == 0) {
			fprintf(stderr, "[%.2f] %i / %i\n", (realtime()-start), counter2, counter);
		}
	}
	kseq_destroy(fqrec1); gzclose(fp1);
	if (file2) {
		kseq_destroy(fqrec2);
		gzclose(fp2);
	} else {
		fqrec2->f = tmp;		// Take care of the original kstream, see above.
		kseq_destroy(fqrec2);   // TODO Check if we need to free or delete tmp
	}
    fprintf(stderr, "[%.2f] Recruited %i reads from %i reads\n", (realtime()-start), counter2, counter);
}

void hash_file(const char *file, const task_t task) {
    if (task == WHITE) {
        fprintf(stderr, "[%.2f] Now whitelisting from %s\n", (realtime()-start), file); // TODO Some fancy foo ? bar : baz
    } else if (task == IGNORE) {
        fprintf(stderr, "[%.2f] Now ignorelisting from %s\n", (realtime()-start), file); // TODO Or move back down
    } else {
        fprintf(stderr, "[%.2f] I don't know what I just did.\n", (realtime()-start)); // TODO assert
        exit(1);
    }

	gzFile fp = strcmp(file, "-") ? gzopen(file, "r") : gzdopen(fileno(stdin), "r");
	kseq_t *seq = kseq_init(fp);

    static int readcount, kmercount;
	while (kseq_read(seq) >= 0) {
        readcount += 1;
		kmercount += hash_read(seq, task);
	}
	kseq_destroy(seq);
	gzclose(fp);

    if (task == WHITE) {
        fprintf(stderr, "[%.2f] Whitelisted %i %i-mers from %i reads\n", (realtime()-start), kmercount, k, readcount);
    } else if (task == IGNORE) {
        fprintf(stderr, "[%.2f] Ignorelisted %i %i-mers from %i reads\n", (realtime()-start), kmercount, k, readcount);
    } else {
        fprintf(stderr, "[%.2f] I don't know what I just did.\n", (realtime()-start));
        exit(1);
    }
}

// TODO -o output file
// TODO tricky Ns with -v
static int usage() {
	fprintf(stderr, "\n%s (%s) version %s by %s (%s)\n\n", PROGRAM, DESCR, VERSION, AUTHOR, CONTACT);

	fprintf(stderr, "Usage: %s [options] <seed> <in1.fq> [in2.fq]\n\n", PROGRAM);

	fprintf(stderr, "       <seed>     recruitment seed (e.g. single-cell reads or reference sequence)\n");
	fprintf(stderr, "       <in1.fq>   read pool to recruit from (e.g. metagenomic reads)\n");
	fprintf(stderr, "       [in2.fq]   optionally: #2 mates, if <in1.fq> contains #1 mates\n\n");

	fprintf(stderr, "Options:\n\n");

	fprintf(stderr, "       -p         <in1.fq> consists of interleaved mates (i.e. shuffled reads)\n");
	fprintf(stderr, "       -x         enable on-the-fly recruitment seed expansion (GREEDY and EXPERIMENTAL)\n\n");

	fprintf(stderr, "       -k INT     k-mer size (15 <= k <= 21), requires 4^(k-1) bytes of RAM (15:256M, 16:1G, 17:4G, 18:16G, ...) [%i]\n\n", k);

	fprintf(stderr, "       -e FLOAT   default: automatically determine the expected error rate based on Phred quality scores with an upper limit [%.2f]\n", max_error);
	fprintf(stderr, "       -d INT     alternatively: min. allowed edit distance (q-gram lemma applied: hits_per_read = readlength - k+1 - d*k) [%i]\n", 0); // TODO not really edit distance, but smth. else
	fprintf(stderr, "       -n INT     alternatively: fixed number of required hits per read [%i]\n\n", 0);

	fprintf(stderr, "       -i FILE    ignore sequence patterns in FILE (e.g. adapter sequences) [null]\n");
    fprintf(stderr, "       --strict   require hits in both mates (?)\n"); // TODO
	fprintf(stderr, "       -c         print only a count of recruited sequences\n");
	fprintf(stderr, "       -v         select non-matching sequences\n\n");

	fprintf(stderr, "       -h         minimal help message\n");
	fprintf(stderr, "       --help     extended help message\n\n"); // TODO 

	return 42;
}

int main(int argc, char *argv[]) {
	start = realtime();
	char *seed = 0, *in1 = 0, *in2 = 0;
	char *ignore = 0;
	int c;
	// TODO getopt_long ?
	while((c = getopt(argc, argv, "pxk:e:d:n:i:cvh")) != -1) {
		switch (c) {
			case 'p':
				flag |= FLAG_P;
				break;
			case 'x':
				flag |= FLAG_X;
				break;
			case 'k':
				k = atoi(optarg);
				if (k < 15) k = 15;
				if (k > 21) k = 21; // Yes, that's 1TB of RAM
				break;
			case 'e':
				max_error = atof(optarg);
				if (max_error > 1.0) max_error = 1.0; // TODO Print warning if max_error > some threshold
				break;
			case 'd':
				if (flag&FLAG_N) return usage();
				flag |= FLAG_D;
				d = atoi(optarg);
				if (d < 0) d = 0;
				break;
			case 'n':
				if (flag&FLAG_D) return usage();
				flag |= FLAG_N;
				d = atoi(optarg);
				if (d < 0) d = 0;
				break;
			case 'i':
				flag |= FLAG_I;
				ignore = optarg;
				break;
			case 'c':
				flag |= FLAG_C;
				break;
			case 'v':
				flag |= FLAG_V;
				break;
			case 'h': // TODO zero exit code on usage explicitely called
			default:
				return usage();
		}
	}
    // TODO This whole flag stuff is rather ugly...

	// TODO Re-write more elegant
	if (optind + 2 == argc) {
		seed = argv[optind];
		in1 = argv[optind+1];
	} else if (optind + 3 == argc && !(flag&FLAG_P)) { // TODO Rather print warning and ignore in2?
		seed = argv[optind];
		in1 = argv[optind+1];
		in2 = argv[optind+2];
	} else return usage();

	fprintf(stderr, "[%.2f] %s (%s) version %s\n", (realtime()-start), PROGRAM, DESCR, VERSION);
	fprintf(stderr, "[%.2f] Executing '%s ", (realtime()-start), PROGRAM);
	if (flag&FLAG_P) fprintf(stderr, "-p ");
	if (flag&FLAG_X) fprintf(stderr, "-x ");
	fprintf(stderr, "-k %i ", k);
	if (flag&FLAG_D) fprintf(stderr, "-d %i ", d);
	if (flag&FLAG_N) fprintf(stderr, "-n %i ", d);
	if (flag&FLAG_I) fprintf(stderr, "-d %s ", ignore);
	if (flag&FLAG_C) fprintf(stderr, "-c ");
	if (flag&FLAG_V) fprintf(stderr, "-v ");
	fprintf(stderr, "%s %s", seed, in1);
	if (in2) fprintf(stderr, " %s", in2);
	fprintf(stderr, "'\n");

	fprintf(stderr, "[%.2f] Allocating %llu MB of memory to store k-mer occurrences\n", (realtime()-start), (1ULL<<(k<<1)>>22));
	kmers = (uint64_t*) calloc(1ULL<<(k<<1)>>5, sizeof(uint64_t));
	if (kmers == NULL) {
		fprintf(stderr, "[%.2f] Failed to allocate enough memory!\n", (realtime()-start));
		return 1;
	}
	if (ignore) {
        hash_file(ignore, IGNORE);
    }
	if (seed) {
        hash_file(seed, WHITE);
    }
	if (in1 && in2) {
		recruit_pe(in1, in2, RECRUIT);
    }
	else if (flag&FLAG_P) {
		recruit_pe(in1, 0, RECRUIT);
    }
	else {
		//TODO recruit_se(in1, RECRUIT);
    }
	free(kmers);
	fprintf(stderr, "[%.2f] Done. Thank you!\n", (realtime()-start));

	return 0;
}
