# Makefile for Bonus Challenges
# This Makefile helps compile and run the bonus challenge simulations

# Basilisk compiler (adjust path if needed)
QCC = qcc

# Compiler flags
CFLAGS = -O2 -Wall

# List of bonus challenges
BONUSES = Bonus1_VaryRe Bonus2_CircularLeaf Bonus3_VaryConcentration \
          Bonus4_InclinedLeaf Bonus5_MultipleLeaves

# Default target: show help
.PHONY: help
help:
	@echo "Bonus Challenges - Makefile"
	@echo "============================"
	@echo ""
	@echo "Available targets:"
	@echo "  make Bonus1_VaryRe.tst           - Run Bonus 1 (Vary Reynolds Number)"
	@echo "  make Bonus2_CircularLeaf.tst     - Run Bonus 2 (Circular Leaf)"
	@echo "  make Bonus3_VaryConcentration.tst - Run Bonus 3 (Vary Concentration)"
	@echo "  make Bonus4_InclinedLeaf.tst     - Run Bonus 4 (Inclined Leaf)"
	@echo "  make Bonus5_MultipleLeaves.tst   - Run Bonus 5 (Multiple Leaves)"
	@echo "  make all                          - Run all bonus challenges"
	@echo "  make clean                        - Clean all generated files"
	@echo ""
	@echo "Example usage:"
	@echo "  make Bonus1_VaryRe.tst"
	@echo ""

# Run all bonus challenges
.PHONY: all
all: $(addsuffix .tst, $(BONUSES))
	@echo "All bonus challenges completed!"

# Generic rule to compile and run any bonus challenge
%.tst: %.c
	@echo "Compiling and running $<..."
	@mkdir -p $*
	@$(QCC) $(CFLAGS) -o $*/$* $<
	@cd $* && ./$* 2>&1 | tee run.log
	@echo "Done! Results in $*/ folder"
	@echo ""

# Clean all generated files
.PHONY: clean
clean:
	@echo "Cleaning bonus challenge outputs..."
	@rm -rf $(BONUSES)
	@echo "Clean complete!"

# Individual shortcuts (for convenience)
.PHONY: bonus1 bonus2 bonus3 bonus4 bonus5
bonus1: Bonus1_VaryRe.tst
bonus2: Bonus2_CircularLeaf.tst
bonus3: Bonus3_VaryConcentration.tst
bonus4: Bonus4_InclinedLeaf.tst
bonus5: Bonus5_MultipleLeaves.tst
