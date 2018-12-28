#= Makefile for deepconv/arare source code in setup.
#
# Authors::   ODAKA Masatsugu, SUGIYAMA Ko-ichiro
# Version::   $Id: Makefile,v 1.10 2014/05/28 15:31:50 sugiyama Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:../../COPYRIGHT]
#
#
#== History
#
# * 2011/06/20 (TAKAHASHI YOSHIYUKI): Add constants0.o
# * 2008/06/19 (ODAKA Masatsugu): Update
# * 2007/08/07 (ODAKA Masatsugu): Update
# * 2006/09/28 (ODAKA Masatsugu): Update
# * 2006/09/12 (ODAKA Masatsugu): Update
# * 2005/01/31 (ODAKA Masatsugu): Update
# * 2005/01/28 (SUGIYAMA Ko-ichiro):
#
########################################################################

include ../../Config.mk

2DOBJS = \
	mpi_wrapper.o \
	argset.o \
	namelist_util.o \
	gridset_2d.o \
	timeset.o \
	clockset.o \
	axesset.o \
	constants0.o \
	constants.o \
	composition.o \
	fileset.o \
	basicset_2d.o 

OBJS = \
	mpi_wrapper.o \
	argset.o \
	namelist_util.o \
	gridset.o \
	timeset.o \
	clockset.o \
	axesset.o \
	constants0.o \
	constants.o \
	composition.o \
	fileset.o \
	basicset.o 

MAIN = test-setup.o

all: $(OBJS)  
	$(AR) $(ARFLAGS) $(LIBDIR)/lib$(LIBNAME).a ./*.o ; \
	$(CP) ./*.mod  $(MODDIR); \

arare4_3d: $(OBJS)
	$(AR) $(ARFLAGS) $(LIBDIR)/lib$(LIBNAME).a ./*.o ; \
	$(CP) ./*.mod  $(MODDIR); \

arare4_2d: $(2DOBJS)  
	$(AR) $(ARFLAGS) $(LIBDIR)/lib$(LIBNAME).a ./*.o ; \
	$(CP) ./*.mod  $(MODDIR); \

test: all $(MAIN)
	$(LINKF) $(MAIN) -L$(LIBDIR) -l$(LIBNAME) -I$(MODDIR) $(SYSLDFLAGS) $(SYSLDLIBS)

%.o: %.F90
	$(FC) $(SYSFFLAGS) $(CPPFLAGS) -c $< 

%.o: %.f90
	$(FC) $(SYSFFLAGS) -c $< 

clean: 
	-$(RM) *.o *.mod  

clean.test: clean
	-$(RM) test *~  
