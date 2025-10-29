CFLAGS += -O2
include $(BASILISK)/Makefile.defs

# Clean up simulation output files
.PHONY: cleanout
cleanout:
	rm -rf *.mp4 dump-00* *.tags.d *.s.d W12 a.out *.dims

# Clean everything (compiled files + output files)
.PHONY: cleanall
cleanall: clean cleanout
