/* The MIT License

   Copyright (c) 2012-2013 Andreas Bremges <abremges@cebitec.uni-bielefeld.de>

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
#define VERSION "0.6.1"
#define AUTHOR  "Andreas Bremges"
#define CONTACT "andreas@cebitec.uni-bielefeld.de"

#include <stdint.h> // uint64_t
#include <stdio.h> // printf
#include <sys/time.h> // gettimeofday (obsolescent)
#include <unistd.h> // getopt
#include <zlib.h> // gzip

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

#define LOGMOD 100000

static double start;
static uint64_t *kmers;
static int k = 16;
static int d = 1;
static int flag = 0;

// A = 00, T = 01, C = 10, G = 11, i.e. XOR with 1 -> compl.
unsigned char seq_fwd_table[128] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

double realtime() {
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

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

int hash_read(const kseq_t *seq, const task_t task) {
	int index = 0;
	int khits = 0;
	
	uint64_t forward = 0;
	uint64_t reverse = 0;
	
	int i;
	for (i = 0; i < seq->seq.l; ++i) {
		index++;
		int c = seq_fwd_table[(int) seq->seq.s[i]];
		if (c < 4) {
			uint64_t mask = (1ULL<<(k<<1))-1;
			forward = (((forward << 2) | c) & mask); // see below
			//if (mask > 0) forward &= mask; <- only a problem if k = 32
			uint64_t shft = (k<<1)-2;
			reverse = ((reverse >> 2) | (((uint64_t) (c^1)) << shft));
		} else {
			index = 0; //return 0;
		}
		if (index >= k) {
			if (task > 0) { // fill hash
				uint64_t z = ((forward < reverse) ? forward : reverse);
				if ((kmers[z>>5]>>((z&31)<<1)&3) < task) {
					kmers[z>>5] |= (task << ((z&31)<<1));
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
void interleaved(const char *file1, const char *file2, const task_t task) {
	
	gzFile fp1 = strcmp(file1, "-") ? gzopen(file1, "r") : gzdopen(fileno(stdin), "r");
	kseq_t *fqrec1 = kseq_init(fp1);
	
	gzFile fp2 = NULL;
	kseq_t *fqrec2 = NULL;
	kstream_t *tmp = NULL;
	
	if (file2) {
		fp2 = gzopen(file2, "r");
		fqrec2 = kseq_init(fp2);
	} else {
		fqrec2 = kseq_init(fp1);
		tmp = fqrec2->f;		// Maybe this can be handled more elegantly?
		fqrec2->f = fqrec1->f;	// fqrec1 and fqrec2 now share one kstream.
	}
	
	static int counter, counter2;
	while (kseq_read(fqrec1) >= 0 && kseq_read(fqrec2) >= 0) {
		counter += 2;
		int n1, n2;
		if ((n1 = fqrec1->seq.l - k+1 - d*k) > 0 && (n2 = fqrec2->seq.l - k+1 - d*k) > 0) {
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
		kseq_destroy(fqrec2);
	}
	if (counter%LOGMOD != 0) {
		fprintf(stderr, "[%.2f] %i / %i\n", (realtime()-start), counter2, counter);
	}
}

void recruit_se(const char *file1, const task_t task) {
	
	gzFile fp1 = strcmp(file1, "-") ? gzopen(file1, "r") : gzdopen(fileno(stdin), "r");
	kseq_t *fqrec1 = kseq_init(fp1);
	
	static int counter, counter2;
	while (kseq_read(fqrec1) >= 0) {
		counter += 1;
		int n1;
		if ((n1 = fqrec1->seq.l - k+1 - d*k) > 0) {
			int khits1;
			if ((khits1 = hash_read(fqrec1, task)) >= n1) {
				if (!(flag&FLAG_V)) {
					counter2 += 1;
					stk_printseq(fqrec1);
				}
				if (flag&FLAG_X) {
					hash_read(fqrec1, WHITE);
				}
			} else if (flag&FLAG_V) {
				counter2 += 1;
				stk_printseq(fqrec1);
			}
		}
		if (counter%LOGMOD == 0) {
			fprintf(stderr, "[%.2f] %i / %i\n", (realtime()-start), counter2, counter);
		}
	}
	kseq_destroy(fqrec1); gzclose(fp1);
	if (counter%LOGMOD != 0) {
		fprintf(stderr, "[%.2f] %i / %i\n", (realtime()-start), counter2, counter);
	}
}

void hash_file(const char *file, const task_t task) {
	gzFile fp = strcmp(file, "-") ? gzopen(file, "r") : gzdopen(fileno(stdin), "r");
	kseq_t *seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		hash_read(seq, task);
	}
	kseq_destroy(seq);
	gzclose(fp);
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
	
	fprintf(stderr, "       -k INT     k-mer size (15 <= k <= 18), requires 4^(k-1) bytes of RAM (15:256M, 16:1G, 17:4G, 18:16G) [%i]\n\n", k);
	
	fprintf(stderr, "       -d INT     min. allowed edit distance (q-gram lemma applied: hits_per_read = readlength - k+1 - d*k) [%i]\n", d);
	fprintf(stderr, "       -n INT     alternatively: fixed number of required hits per read [%i]\n\n", 0);
//	fprintf(stderr, "       -r FLOAT   alternatively: ratio of known k-mers per read [%.2f]\n\n", 0.0);
	
	fprintf(stderr, "       -i FILE    ignore sequence patterns in FILE (e.g. adapter sequences) [null]\n");
	fprintf(stderr, "       -c         print only a count of recruited sequences\n");
	fprintf(stderr, "       -v         select non-matching sequences\n\n");
	
	return 42;
}

int main(int argc, char *argv[]) {
	start = realtime();
	char *seed = 0, *in1 = 0, *in2 = 0;
	char *ignore = 0;
	int c;
	while((c = getopt(argc, argv, "pxk:d:n:i:cv")) != -1) {
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
				if (k > 18) k = 18;
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
			default:
				return usage();
		}
	}
	if (!(flag&FLAG_N)) flag |= FLAG_D;

	// TODO Re-write more elegant
	if (optind + 2 == argc) {
		seed = argv[optind];
		in1 = argv[optind+1];
	} else if (optind + 3 == argc && !(flag&FLAG_P)) {
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
		fprintf(stderr, "[%.2f] Now ignorelisting from %s\n", (realtime()-start), ignore);
		hash_file(ignore, IGNORE);
		fprintf(stderr, "[%.2f] Ignorelisted xx %i-mers from yy reads\n", (realtime()-start), k);
	}
	if (seed) {
		fprintf(stderr, "[%.2f] Now whitelisting from %s\n", (realtime()-start), seed);
		hash_file(seed, WHITE);
		fprintf(stderr, "[%.2f] Whitelisted xx %i-mers from yy reads\n", (realtime()-start), k);
	}
	if (in1 && in2) {
		fprintf(stderr, "[%.2f] Now recruiting pairs from %s and %s\n", (realtime()-start), in1, in2);
		interleaved(in1, in2, RECRUIT);
		fprintf(stderr, "[%.2f] Recruited xx reads from yy reads\n", (realtime()-start));
	} else if (flag&FLAG_P) {
		fprintf(stderr, "[%.2f] Now recruiting pairs from %s\n", (realtime()-start), in1);
		interleaved(in1, 0, RECRUIT);
		fprintf(stderr, "[%.2f] Recruited xx reads from yy reads\n", (realtime()-start));
	} else {
		fprintf(stderr, "[%.2f] Now recruiting from %s\n", (realtime()-start), in1);
		recruit_se(in1, RECRUIT);
		fprintf(stderr, "[%.2f] Recruited xx reads from yy reads\n", (realtime()-start));
	}
	free(kmers);
	fprintf(stderr, "[%.2f] Done. Thank you!\n", (realtime()-start));

	return 0;
}
