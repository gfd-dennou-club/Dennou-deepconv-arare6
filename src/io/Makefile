#= Makefile for deepconv/arare source code in io.
#
# Authors::   ODAKA Masatsugu, SUGIYAMA Ko-ichiro
# Version::   $Id: Makefile,v 1.11 2010-08-11 07:34:19 sugiyama Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:../../COPYRIGHT]
#
#
#== History
#
# * 2008/06/19 (ODAKA Masatsugu): Update
# * 2007/08/03 (ODAKA Masatsugu): Update
# * 2006/09/28 (ODAKA Masatsugu): Update
# * 2006/09/12 (ODAKA Masatsugu): Update
# * 2005/04/16 (SUGIYAMA Ko-ichiro):
#
########################################################################

include ../../Config.mk

.SUFFIXES: .f90 .o 

2DOBJS = \
	historyfileio_2d.o \
	restartfileio_2d.o \
	basicfileio_2d.o

OBJS = \
	historyfileio.o \
	restartfileio.o \
	restartfileio_ih1998.o \
	basicfileio.o \
	arare4restartfileio.o \
	arare4initfileio.o 

all: $(OBJS) 
	$(AR) $(ARFLAGS) $(LIBDIR)/lib$(LIBNAME).a ./*.o ; \
	$(CP) ./*.mod  $(MODDIR); \


.f90.o: 
	$(FC) $(SYSFFLAGS) -c $< 

clean: 
	-$(RM) *.o *.mod

