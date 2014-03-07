CC=gcc
CFLAGS=-Wall -O3

kgrep:kgrep.c kseq.h
	$(CC) $(CFLAGS) kgrep.c -o $@ -lz

clean:
	rm -fr gmon.out *.o ext/*.o a.out kgrep *~ *.a *.dSYM session*
