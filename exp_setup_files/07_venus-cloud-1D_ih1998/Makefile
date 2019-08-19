#= Makefile for deepconv/arare source code in chemdata.
#
# Authors::   ODAKA Masatsugu, SUGIYAMA Ko-ichiro
# Version::   $Id: Makefile,v 1.10 2009-03-05 04:37:21 sugiyama Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:../../COPYRIGHT]
#
#
#== History
#
# * 2006/09/28 (ODAKA Masatsugu): Update
# * 2006/09/12 (ODAKA Masatsugu): Update
# * 2006/01/09 (SUGIYAMA Ko-ichiro):
#
########################################################################

include ../../Config.mk

.SUFFIXES : .f90 .o

MAIN = \
	main.o

all:   $(MAIN)
	$(LINKF) $(OBJS) $(MAIN) -o main \
	-L$(LIBDIR) -l$(LIBNAME) -I$(MODDIR)

clean: 
	-$(RM) *.o main

.f90.o: 
	$(FC) $(SYSFFLAGS) -c $< 

