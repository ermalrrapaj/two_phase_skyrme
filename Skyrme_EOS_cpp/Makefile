CC=gcc
CFLAGS=-g -O3 -fopenmp 
CXX=g++
CPPFLAGS=-g -O3 -std=c++11 -fopenmp 
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
#INCLUDE += $(HDF5INCS)  
#LIBS += $(HDF5LIBS)

LIBS += -lgsl -lgslcblas -lm 

# Source files
UTILDIR=./
UTILS:= 
UTILS:=$(addprefix $(UTILDIR)/,$(UTILS)) 

EOSDIR=./
EOSSRCS:= funct_fermi1.f90 EOSData.cpp EOSSkyrme.cpp GibbsPhaseConstruct.cpp EOSSingleNucleus.cpp
EOSSRCS:=$(addprefix $(EOSDIR)/,$(EOSSRCS)) 

SRCS_ := $(EOSSRCS) $(UTILS) 
OBJS_ := ${SRCS_:.F90=.o}
OBJS_ := ${OBJS_:.f90=.o}
OBJS_ := ${OBJS_:.cpp=.o}
OBJS = $(patsubst %,$(OBJDIR)/%,$(OBJS_)) 
SRCS = $(patsubst %,$(SRCDIR)/%,$(SRCS_)) 
 
all: create_dirs skyrme gibbs_test

skyrme: $(OBJS) $(OBJDIR)/skyrme_exe.o
	$(CXX) $(INCLUDE) -o $(EXEDIR)/$@ $^ $(LIBS)

gibbs_test: $(OBJS) $(OBJDIR)/gibbs_test.o
	$(CXX) $(INCLUDE) -o $(EXEDIR)/$@ $^ $(LIBS)

create_dirs:
	@ if test ! -d $(OBJDIR); then mkdir $(OBJDIR); fi
	@ if test ! -d $(LIBDIR); then mkdir $(LIBDIR); fi
	@ if test ! -d $(EXEDIR); then mkdir $(EXEDIR); fi
	@ if test ! -d $(INCDIR); then mkdir $(INCDIR); fi
	@ if test ! -d $(OBJDIR)/$(NRDIR); then mkdir $(OBJDIR)/$(NRDIR); fi
	@ if test ! -d $(OBJDIR)/$(EOSDIR); then mkdir $(OBJDIR)/$(EOSDIR); fi
	@ if test ! -d $(OBJDIR)/$(UTILDIR); then mkdir $(OBJDIR)/$(UTILDIR); fi

clean:
	find $(OBJDIR) -name \*.[o,a] -delete	
	find $(INCDIR) -name \*.mod -delete	

$(OBJDIR)/$(notdir %.o): $(SRCDIR)/%.c
	@echo Compiling $< to $(OBJDIR)/$(notdir $@) 
	@$(CC) -c $(CFLAGS) $(INCLUDE) -o $@ $< $(LIBS)

$(OBJDIR)/$(notdir %.o): $(SRCDIR)/%.cpp
	@echo Compiling $< to $(OBJDIR)/$(notdir $@) 
	@$(CXX) -c $(CPPFLAGS) $(INCLUDE) -o $@ $< $(LIBS)

$(OBJDIR)/$(notdir %.o): $(SRCDIR)/%.f90 
	@echo Compiling $< to $(OBJDIR)/$(notdir $@) 
	@$(F90) -J$(INCDIR) -c $(F90FLAGS) $(INCLUDE) -o $@ $< $(LIBS)

$(OBJDIR)/$(notdir %.o): $(SRCDIR)/%.F90 
	@echo Compiling $< to $@
	@$(F90) -J$(INCDIR) -c $(F90FLAGS) $(INCLUDE) -o $@ $< $(LIBS)
