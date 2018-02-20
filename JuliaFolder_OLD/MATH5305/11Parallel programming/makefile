DEBUG_FLAGS = -fcheck=all
OPT_FLAGS = -O2
FFLAGS = $(OPT_FLAGS)

maxeig: maxeig.f08
	gfortran $(FFLAGS) -o $@ $<

clean:
	rm -f maxeig
