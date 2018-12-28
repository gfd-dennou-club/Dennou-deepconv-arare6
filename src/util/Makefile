#= Makefile for deepconv/arare source code in util.
#
# Authors::   ODAKA Masatsugu, SUGIYAMA Ko-ichiro
# Version::   $Id: Makefile,v 1.11 2009-03-05 04:37:24 sugiyama Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:../../COPYRIGHT]
#
#
#== History
#
# * 2006/09/28 (ODAKA Masatsugu): Update
# * 2006/09/12 (ODAKA Masatsugu): Update
# * 2005/04/16 (SUGIYAMA Ko-ichiro):
#
########################################################################

include ../../Config.mk

.SUFFIXES: .f90 .o

2DOBJS = \
	average_2d.o \
	differentiate_center2_2d.o \
	differentiate_center4_2d.o \
	boundary.o \
	damping_2d.o \
	cflcheck.o \
	fillnegative_2d.o \
	numdiffusion4th_2d.o 

3DOBJS = \
	average.o \
	differentiate_center2.o \
	differentiate_center4.o \
	setmargin.o \
	damping.o \
	cflcheck.o \
	fillnegative.o \
	numdiffusion4th.o 

OBJS = \
	average.o \
	differentiate_center2.o \
	differentiate_center4.o \
	setmargin.o \
	damping.o \
	cflcheck.o \
	fillnegative.o \
	energymonit.o \
	dexnerdt.o

all: $(OBJS) 
	$(AR) $(ARFLAGS) $(LIBDIR)/lib$(LIBNAME).a ./*.o ; \
	$(CP) ./*.mod  $(MODDIR); \

arare4_2d: $(2DOBJS) 
	$(AR) $(ARFLAGS) $(LIBDIR)/lib$(LIBNAME).a ./*.o ; \
	$(CP) ./*.mod  $(MODDIR); \

arare4_3d: $(3DOBJS) 
	$(AR) $(ARFLAGS) $(LIBDIR)/lib$(LIBNAME).a ./*.o ; \
	$(CP) ./*.mod  $(MODDIR); \

.f90.o: 
	$(FC) $(SYSFFLAGS) -c $< 


clean: 
	-$(RM) *.o *.mod 
