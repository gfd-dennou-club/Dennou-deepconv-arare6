#= Makefile for deepconv/arare source code in physics.
#
# Authors::   ODAKA Masatsugu, SUGIYAMA Ko-ichiro
# Version::   $Id: Makefile,v 1.12 2009-03-05 04:37:23 sugiyama Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:../../COPYRIGHT]
#
#
#== History
#
# * 2006/09/28 (ODAKA Masatsugu): Update
# * 2006/09/12 (ODAKA Masatsugu): Update
# * 2005/01/31 (ODAKA Masatsugu): Update
# * 2005/01/28 (SUGIYAMA Ko-ichiro):
#
########################################################################

include ../../Config.mk

.SUFFIXES: .f90 .o 

OBJS = \
	surfaceflux_baker1998.o \
	surfaceflux_diff.o \
	surfaceflux_bulk.o \
	surfaceflux_const.o 

all: $(OBJS) 
	$(AR) $(ARFLAGS) $(LIBDIR)/lib$(LIBNAME).a ./*.o ; \
	$(CP) ./*.mod  $(MODDIR); \

.f90.o: 
	$(FC) $(SYSFFLAGS) -c $< 


clean: 
	-$(RM) *.o *.mod 


clean.test: clean
	-$(RM) test *~
