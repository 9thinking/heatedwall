# Makefile for Basilisk droplet impact simulation

# Basilisk installation directory (modify if needed)
BASILISK = $(HOME)/basilisk/src

# Compiler settings
CC = gcc -std=c99
CFLAGS = -O2 -Wall -I$(BASILISK)
LDFLAGS = -lm

# Basilisk preprocessor
QCC = qcc -autolink -O2 -Wall

# Targets
all: droplet_impact

droplet_impact: droplet_impact.c
	$(QCC) droplet_impact.c -o droplet_impact -lm

# Clean build files
clean:
	rm -f droplet_impact *.o *.s
	rm -f *.ppm *.gif *.mp4
	rm -f *.dat
	rm -f _*.c

# Clean all output files
cleanall: clean
	rm -f interface-*.dat field-*.dat
	rm -f quantitative_output.dat
	rm -f log

# Run simulation with default parameters
run: droplet_impact
	./droplet_impact

# Run with higher resolution
run-hr: droplet_impact
	./droplet_impact 10 12

# Run with custom impact velocity (example: 2 m/s)
run-fast: droplet_impact
	./droplet_impact 9 11 2.0

# Generate plots from output data
plot:
	@echo "Generating plots with gnuplot..."
	gnuplot plot_results.gp

# Help message
help:
	@echo "Makefile for Basilisk droplet impact simulation"
	@echo ""
	@echo "Targets:"
	@echo "  make all       - Compile the simulation"
	@echo "  make run       - Run with default parameters"
	@echo "  make run-hr    - Run with higher resolution"
	@echo "  make run-fast  - Run with 2 m/s impact velocity"
	@echo "  make plot      - Generate plots from output data"
	@echo "  make clean     - Remove compiled files"
	@echo "  make cleanall  - Remove all output files"
	@echo ""
	@echo "Command line arguments:"
	@echo "  ./droplet_impact [LEVEL] [MAXLEVEL] [U0]"
	@echo "    LEVEL    - Base grid level (default: 9)"
	@echo "    MAXLEVEL - Maximum refinement level (default: 11)"
	@echo "    U0       - Impact velocity in m/s (default: 1.0)"

.PHONY: all clean cleanall run run-hr run-fast plot help
