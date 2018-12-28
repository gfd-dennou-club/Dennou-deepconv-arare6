#= Makefile for deepconv/arare source in main.
#
# Authors::   ODAKA Masatsugu, SUGIYAMA Ko-ichiro
# Version::   $Id: Makefile,v 1.20 2011-02-28 12:00:24 sugiyama Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:../../COPYRIGHT]
#
#
#== History
#
# * 2018/12/28 (SUGIYAMA Ko-ichiro): Update
# * 2006/09/28 (ODAKA Masatsugu): Update
# * 2006/09/12 (ODAKA Masatsugu): Update
# * 2005/04/10 (SUGIYAMA Ko-ichiro):
#
########################################################################

include ../../Config.mk


.SUFFIXES: .f90 .o 

OBJS = arare.o arare_init-data.o

all: arare init

arare: $(OBJS) $(LIBDIR)/lib$(LIBNAME).a 
	$(LINKF) arare.o -o arare \
	-L$(LIBDIR) -l$(LIBNAME) -I$(MODDIR) \
	$(SYSLDFLAGS) $(SYSLDLIBS)

	$(CP) arare $(BINDIR)

init: $(OBJS) $(LIBDIR)/lib$(LIBNAME).a 
	$(LINKF) arare_init-data.o -o arare_init-data \
	-L$(LIBDIR) -l$(LIBNAME) -I$(MODDIR) \
	$(SYSLDFLAGS) $(SYSLDLIBS)

	$(CP) arare_init-data $(BINDIR)

.f90.o: 
	$(FC) $(SYSFFLAGS) -c $< 


clean:
	-$(RM) -f *.o arare arare_init-data *~

