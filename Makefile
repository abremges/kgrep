CC=gcc
CFLAGS=-Wall -pedantic -std=gnu99 -O3

kgrep:kgrep.c kseq.h
	$(CC) $(CFLAGS) kgrep.c -o $@ -lz -lm

clean:
	rm -fr gmon.out *.o ext/*.o a.out kgrep *~ *.a *.dSYM session*
