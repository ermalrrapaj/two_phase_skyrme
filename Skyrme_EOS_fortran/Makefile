F90=gfortran
F90FLAGS=-g -O3 -fopenmp 
INCLUDE = 
LIBS = 

# Local directory structure
OBJDIR = obj
SRCDIR = src
LIBDIR = lib
INCDIR = include
EXEDIR = exe
INCLUDE += -I$(INCDIR) 

# HDF5 libraries
HDF5INCS=-I/usr/include
HDF5LIBS=-L/usr/lib/x86_64-linux-gnu/ -lhdf5 -lhdf5_fortran -lhdf5 -lz -lhdf5_hl -lhdf5hl_fortran
INCLUDE += $(HDF5INCS)  
LIBS += $(HDF5LIBS)

# Source files
NRDIR=nr_src
NRSRCS := nrtype.f90 nrutil.f90 nr_interface.f90 nr.f90 fdjac.f90 fmin.f90 lnsrch.f90 ludcmp.f90 lubksb.f90 newt.f90
NRSRCS:=$(addprefix $(NRDIR)/,$(NRSRCS)) 

UTILDIR=util
UTILS:= h5_output.F90 bisect.f90 funct_fermi1.f90 
UTILS:=$(addprefix $(UTILDIR)/,$(UTILS)) 

EOSDIR=eos
EOSSRCS:= eos_com.F90 eos_skyrme.F90 eos_nse.F90 eos_gibbs_region.F90 LS_construction.F90 Tabular_EOS.F90 
EOSSRCS:=$(addprefix $(EOSDIR)/,$(EOSSRCS)) 

SRCS_ := $(NRSRCS) $(UTILS) $(EOSSRCS)
OBJS_ := ${SRCS_:.F90=.o}
OBJS_ := ${OBJS_:.f90=.o}
OBJS = $(patsubst %,$(OBJDIR)/%,$(OBJS_)) 
SRCS = $(patsubst %,$(SRCDIR)/%,$(SRCS_)) 
 
all: create_dirs bulk LS gibbs_bound nucleus_properties build_table nse_test

create_dirs:
	if test ! -d $(OBJDIR); then mkdir $(OBJDIR); fi
	if test ! -d $(LIBDIR); then mkdir $(LIBDIR); fi
	if test ! -d $(EXEDIR); then mkdir $(EXEDIR); fi
	if test ! -d $(INCDIR); then mkdir $(INCDIR); fi
	if test ! -d $(OBJDIR)/$(NRDIR); then mkdir $(OBJDIR)/$(NRDIR); fi
	if test ! -d $(OBJDIR)/$(EOSDIR); then mkdir $(OBJDIR)/$(EOSDIR); fi
	if test ! -d $(OBJDIR)/$(UTILDIR); then mkdir $(OBJDIR)/$(UTILDIR); fi

nse_test: $(OBJS) $(OBJDIR)/nse_test.o
	$(F90) $(F90FLAGS) $(INCLUDE) -o $(EXEDIR)/$@ $^ $(LIBS)

nucleus_properties: $(OBJS) $(OBJDIR)/nucleus_properties.o
	$(F90) $(F90FLAGS) $(INCLUDE) -o $(EXEDIR)/$@ $^ $(LIBS)

tabular_test: $(OBJS) $(OBJDIR)/tabular_test.o 
	$(F90) $(F90FLAGS) $(INCLUDE) -o $(EXEDIR)/$@ $^ $(LIBS)
	 
build_table: $(OBJS) $(OBJDIR)/build_table.o 
	$(F90) $(F90FLAGS) $(INCLUDE) -o $(EXEDIR)/$@ $^ $(LIBS)

bulk: $(OBJS) $(OBJDIR)/get_bulk_properties.o
	$(F90) $(F90FLAGS) $(INCLUDE) -o $(EXEDIR)/$@ $^ $(LIBS)

gibbs_bound: $(OBJS) $(OBJDIR)/get_gibbs_boundary.o
	$(F90) $(F90FLAGS) $(INCLUDE) -o $(EXEDIR)/$@ $^ $(LIBS)

LS: $(OBJS) $(OBJDIR)/LSexe.o
	$(F90) $(F90FLAGS) $(INCLUDE) -o $(EXEDIR)/$@ $^ $(LIBS)

h5_writer_test: $(OBJS) $(OBJDIR)/h5_writer_test.o
	$(F90) $(F90FLAGS) $(INCLUDE) -o $(EXEDIR)/$@ $^ $(LIBS)

skyrme: $(OBJS) $(OBJDIR)/Skyrme_EOS.o $(OBJDIR)/skyrme_exe.o
	$(F90) $(F90FLAGS) $(INCLUDE) -o $(EXEDIR)/$@ $^ $(LIBS)

clean:
	find $(OBJDIR) -name \*.[o,a] -delete	
	find $(INCDIR) -name \*.mod -delete	

$(OBJDIR)/$(notdir %.o): $(SRCDIR)/%.f90 
	@echo Compiling $< to $(OBJDIR)/$(notdir $@) 
	@$(F90) -J$(INCDIR) -c $(F90FLAGS) $(INCLUDE) -o $@ $< $(LIBS)

$(OBJDIR)/$(notdir %.o): $(SRCDIR)/%.F90 
	@echo Compiling $< to $@
	@$(F90) -J$(INCDIR) -c $(F90FLAGS) $(INCLUDE) -o $@ $< $(LIBS)
