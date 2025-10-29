CFLAGS += -O2
include $(BASILISK)/Makefile.defs

# Green2D MPI target
Green2D.tst: Green2D.c
	CC99='mpicc -std=c99' qcc -Wall -O2 -disable-dimensions -D_MPI=1 Green2D.c -lm -lglutils -L$$BASILISK/gl -lfb_tiny
	mpiexec -np 4 ./a.out

# Clean up simulation output files
.PHONY: cleanout
cleanout:
	rm -rf *.mp4 dump-00* *.tags.d *.s.d W12 a.out dump-*

# Clean everything (compiled files + output files)
.PHONY: cleanall
cleanall: clean cleanout