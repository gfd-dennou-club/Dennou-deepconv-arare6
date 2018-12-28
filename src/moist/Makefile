#= Makefile for deepconv/arare source code in moist.
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
# * 2008/06/19 (ODAKA Masatsugu): Update
# * 2006/09/28 (ODAKA Masatsugu): Update
# * 2006/09/12 (ODAKA Masatsugu): Update
# * 2005/01/31 (ODAKA Masatsugu): Update
# * 2005/01/28 (SUGIYAMA Ko-ichiro):
#
########################################################################

include ../../Config.mk


.SUFFIXES: .f90 .o 

#OBJS = \
#	chemcalc.o \
#	moistfunc.o    \
#	moistadjust.o \
#	warmrainprm.o \
#	eccm.o        \
#	moistbuoyancy.o  

OBJS = \
	chemcalc.o \
	eccm.o  \
	moistadjust.o \
	cloudphys_k1969.o \
	cloudphys_marscond.o

all: $(OBJS) 
	$(AR) $(ARFLAGS) $(LIBDIR)/lib$(LIBNAME).a ./*.o ; \
	$(CP) ./*.mod  $(MODDIR); \

.f90.o: 
	$(FC) $(SYSFFLAGS) -c $< 

clean: 
	-$(RM) *.o *.mod 


clean.test: clean
	-$(RM) test *~
