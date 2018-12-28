#= Makefile for deepconv/arare source code in env.
#
# Authors::   ODAKA Masatsugu, SUGIYAMA Ko-ichiro
# Version::   $Id: Makefile,v 1.9 2014/07/08 00:59:08 sugiyama Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:../../COPYRIGHT]
#
#
#== History
#
# * 2006/09/28 (ODAKA Masatsugu): Update
# * 2006/09/12 (ODAKA Masatsugu): Update
# * 2006/04/09 (SUGIYAMA Ko-ichiro):
#
########################################################################

include ../../Config.mk


.SUFFIXES: .f90 .o 


OBJS = \
	initialdata_disturb.o\
	initialdata_toon2002.o \
	initialdata_takemi2007.o \
	initialdata_baker1998.o \
	initialdata_sounding.o \
	initialdata_yamasaki1983.o \
	initialdata_skamarock1994.o

all: $(OBJS) 
	$(AR) $(ARFLAGS) $(LIBDIR)/lib$(LIBNAME).a ./*.o ; \
	$(CP) ./*.mod  $(MODDIR); \

.f90.o: 
	$(FC) $(SYSFFLAGS) -c $< 

clean: 
	-$(RM) *.o *.mod
