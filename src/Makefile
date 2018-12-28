#= Makefile for deepconv/arare source code. 
#
# Authors::   ODAKA Masatsugu, SUGIYAMA Ko-ichiro
# Version::   $Id: Makefile,v 1.7 2009-03-05 04:37:21 sugiyama Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:../../COPYRIGHT]
#
#
#== History
#
# * 2006/09/28 (ODAKA Masatsugu): Update
# * 2006/09/12 (ODAKA Masatsugu): Update
# * 2006/03/19 (SUGIYAMA Ko-ichiro): Update
# * 2005/01/27 (ODAKA Masatsugu): Update
# * 2003/02/05 (SUGIYAMA Ko-ichiro):
#
########################################################################

include ../Config.mk

.SUFFIXES: .o .f90

EXEDIR   = main
SHAREDIR = chemdat setup util dynamic moist radiation surface_flux turbulence env io 
SUBDIR   = $(SHAREDIR) $(EXEDIR) 


all: dir
	@for dir in $(SHAREDIR) ; do \
	  cd $$dir ; \
	  $(MAKE) all ; \
	  cd ../ ; \
	done 
	@for dir in $(EXEDIR) ; do \
	  cd $$dir ; \
	  $(MAKE) all ; \
	  cd ../ ; \
	done 

dir: 
	@if ( test ! -d $(BINDIR) ) then  \
	  echo "Making directory $(BINDIR)" ; \
	  mkdir $(BINDIR) ; \
	fi	
	@if (test ! -d $(LIBDIR) ) then  \
	  echo "Making directory $(LIBDIR)" ; \
	  mkdir $(LIBDIR) ; \
	fi
	@if (test ! -d $(MODDIR) ) then  \
	  echo "Making directory $(MODDIR)" ; \
	  mkdir $(MODDIR) ; \
	fi

clean:
	@for dir in $(SUBDIR) ; do \
	  cd $$dir ; \
	  $(MAKE) clean ; \
	  cd ../ ; \
	done 



