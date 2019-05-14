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

OBJS = \
	chemdata.o \
	imamura1998.o

MAIN = \
	test-imamura3.o
#	test-imamura2.o
#	test-imamura.o
#	test-chemdata.o
#	test-zeleznik1999.o

TEST1 = \
	test-imamura-1.o

TEST2 = \
	test-imamura-2.o

TEST3 = \
	test-imamura-3.o

TEST4 = \
	test-imamura-4.o


all: $(OBJS) 
	$(AR) $(ARFLAGS) $(LIBDIR)/lib$(LIBNAME).a ./*.o ; \
	$(CP) ./*.mod  $(MODDIR)

test: $(OBJS) $(MAIN)
	$(LINKF) $(OBJS) $(MAIN) -o test \
	-L$(LIBDIR) -l$(LIBNAME) -I$(MODDIR)

test1: $(OBJS) $(TEST1)
	$(LINKF) $(OBJS) $(TEST1) -o test1 \
	-L$(LIBDIR) -l$(LIBNAME) -I$(MODDIR)

test2: $(OBJS) $(TEST2)
	$(LINKF) $(OBJS) $(TEST2) -o test2 \
	-L$(LIBDIR) -l$(LIBNAME) -I$(MODDIR)

test3: $(OBJS) $(TEST3)
	$(LINKF) $(OBJS) $(TEST3) -o test3 \
	-L$(LIBDIR) -l$(LIBNAME) -I$(MODDIR)

test4: $(OBJS) $(TEST4)
	$(LINKF) $(OBJS) $(TEST4) -o test4 \
	-L$(LIBDIR) -l$(LIBNAME) -I$(MODDIR) 

clean: 
	-$(RM) *.o *.mod 


clean.test: clean 
	-$(RM) test *~

.f90.o: 
	$(FC) $(SYSFFLAGS) -c $< 

