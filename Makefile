# program name
PROG    = ./Geometry


# compiler options
FC= gfortran #-static
FFLAGS = -Og -g -ffree-line-length-none -fbounds-check
LINKER  = $(FC) 

#OPENBLAS=./usr/qcopenblas/
# LIBS     = -I$(OPENBLAS)/include -L$(OPENBLAS)/lib -lopenblas 

# add your own files here:
 SOURCES=\
 geometry.f90\
 main.f90\
 string.f90\
 tools.f90\
 ebe.f
# geometry_prb.f90

SOURCE77=
OF90=$(SOURCES:.f90=.o)
OF77=$(SOURCE77:.f=.o)

OBJS=$(OF90) $(OF77)

%.o: %.f
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@
# link
$(PROG): $(OBJS)
	$(LINKER) $(OBJS) $(LIBS) -o $(PROG)

clean:
	rm -f *.o $(PROG)

tarball:
	tar cvzf Geometry.tgz Makefile $(SOURCES) $(SOURCE77)
