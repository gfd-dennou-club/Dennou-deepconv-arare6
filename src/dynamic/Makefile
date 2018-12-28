#= Makefile for deepconv/arare source code in dynamic.
#
# Authors::   ODAKA Masatsugu, SUGIYAMA Ko-ichiro
# Version::   $Id: Makefile,v 1.15 2009-03-05 04:37:21 sugiyama Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:../../COPYRIGHT]
#
#
#== History
#
# * 2007/08/03 (ODAKA Masatsugu): Update
# * 2006/09/28 (ODAKA Masatsugu): Update
# * 2006/09/12 (ODAKA Masatsugu): Update
# * 2005/01/31 (ODAKA Masatsugu): Update
# * 2005/01/28 (SUGIYAMA Ko-ichiro):
#
########################################################################

include ../../Config.mk


.SUFFIXES: .f90 .o 

3DOBJS = \
	dynfunc_3d.o \
	dynimpfunc_3d.o \
	dynamics_hevi_v1.o \
	dynamics_hevi_v2.o \
	acousticmode_2d.o \
	acousticmode_3d.o \
	acousticmode_std.o \
	advection_center4_2d.o \
	advection_center4_3d.o \
	advection_center4_std.o \
	dynamics_hevi_v3.o 

OBJS = \
	acousticmode_2d.o \
	acousticmode_3d.o \
	acousticmode_std.o \
	advection_center4_2d.o \
	advection_center4_3d.o \
	advection_center4_std.o \
	dynamics_hevi_v3.o 


all: $(OBJS) 
	$(AR) $(ARFLAGS) $(LIBDIR)/lib$(LIBNAME).a ./*.o ; \
	$(CP) ./*.mod  $(MODDIR)

.f90.o: 
	$(FC) $(SYSFFLAGS) -c $< 

dynimpfunc.f90: dynimpfunc.f90.$(LINSOLV)
	ln -s dynimpfunc.f90.$(LINSOLV) $@


dynimpfunc_3d.f90: dynimpfunc_3d.f90.$(LINSOLV)
	ln -s dynimpfunc_3d.f90.$(LINSOLV) $@


clean: 
	-$(RM) *.o *.mod dynimpfunc.f90 dynimpfunc_3d.f90

