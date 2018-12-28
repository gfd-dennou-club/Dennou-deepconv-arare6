#= Makefile for deepconv/arare6 test calculation
#
# Authors::   SUGIYAMA Ko-ichiro
# Version::   $Id: Makefile,v 1.7 2009-03-05 04:37:21 sugiyama Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2015. All rights reserved.
# License::   See COPYRIGHT[link:../../COPYRIGHT]
#
#
#== History
#
# * 2015/07/22 (SUGIYAMA Ko-ichiro): CREATE
#
########################################################################

RUBYDIR = 01_sound-wave 02_advection 02_advection-QMix 03_gravity-wave 04_density-current 05_warm-bubble 06_warm-bubble-kw1978

RUBYEXEC = testrun.rb
RUBY = ruby
RM = rm
XVFB = xvfb-run 

all: 
	@for dir in $(RUBYDIR) ; do \
	  cd $$dir ; \
	  $(XVFB) $(RUBY) $(RUBYEXEC) ; \
          sleep 10; \
	  cd ../ ; \
	done 

clean:
	@for dir in $(RUBYDIR) ; do \
	  cd $$dir ; \
	  $(RM) *conf *png *nc ; \
	  cd ../ ; \
	done 



