# makefile for streamline integrator based on Pollock algorithm
# Daniel W. Meyer, Institute of Fluid Dynamics, ETH Zurich, March 2017
FC = gfortran
OBJS = streamlinesnt3d.o pollock_lib.o streamlinesnt3d_lib.o
CFLAGS = -Og -g -fimplicit-none -Wall -fwhole-file -fopenmp -fcheck=all -fcheck=bounds -ffpe-trap=zero,overflow,invalid
all: streamlinesnt3d
streamlinesnt3d: $(OBJS)
	$(FC) $(OBJS) $(CFLAGS) -o streamlinesnt3d
streamlinesnt3d.o: streamlinesnt3d.f90 pollock_lib.o streamlinesnt3d_lib.o
	$(FC) $(CFLAGS) -c $<
pollock_lib.o: pollock_lib.f90
	$(FC) $(CFLAGS) -c $<
streamlinesnt3d_lib.o: streamlinesnt3d_lib.f90 pollock_lib.f90
	$(FC) $(CFLAGS) -c $<
clean:
	rm -rf $(OBJS) ../streamlinesnt3d.exe pollock_lib.mod streamlinesnt3d_lib.mod streamlinesnt3d.dSYM
