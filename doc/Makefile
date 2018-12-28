#= Makefile for deepconv/arare document source tree.
#
# Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
# Version::   $Id: Makefile,v 1.4 2014/02/21 04:18:12 deepconv Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:../COPYRIGHT]
#
#
#== History
#
# * 2014/02/16 (TAKEHIRO Shin-ichi)
# * 2006/09/29 (ODAKA Masatsugu)
# * 2006/09/16 (ODAKA Masatsugu)
# * 2006/03/19 (SUGIYAMA Ko-ichiro)
# * 2005/01/27 (ODAKA Masatsugu)
#       
########################################################################
#== Settings
#
# [JAPANESE] 各種設定項目

#INCLUDEFILE = ../Config.mk
#include $(INCLUDEFILE)  # Include file
                        # インクルードファイル

# End Settings  [JA] 設定項目ここまで
########################################################################

######################################################################
#== Set subdirectories

SUBDIR  = dai1bu/src dai2bu/src 

######################################################################
#== Rules
#

doc: rd2html rdoc latex2html


rd2html:
	rake
#	$(MAKE) -f Makefile.rd2html
#	cd tutorial ; $(MAKE) -f Makefile.rd2html ; cd ../


rdoc:
	$(MAKE) -f Makefile.rdoc


latex2html:
	@for dir in $(SUBDIR) ; do \
	  cd $$dir ; \
	  $(MAKE) all ; \
	  cd ../.. ; \
	done 

clean: clean.local clean.dir
	-$(RM) *~ *.bak */*~ */*.bak

clean.all: clean.rdoc clean.rd2html clean.latex2html clean.all.dir

clean.local:
	@if [ ! X"$(wildcard *~ *.bak)" = X"" ]; then \
		echo rm $(wildcard *~ *.bak) ;\
		rm $(wildcard *~ *.bak) ;\
	fi

clean.dir:
	@for i in $(SUBDIR) ; do \
	  cd $$i || exit 1; \
	  echo "Cleaning $$i/" ; \
	  $(MAKE) clean || exit 1 ; \
	  cd ../.. ; \
	done

clean.all.dir:
	@for i in $(SUBDIR) ; do \
	  cd $$i || exit 1; \
	  echo "Cleaning $$i/" ; \
	  $(MAKE) clean.all || exit 1 ; \
	  cd ../.. ; \
	done

clean.rd2html:
	rake distclean
#	$(MAKE) clean.all -f Makefile.rd2html
#	cd tutorial ; $(MAKE) clean.all -f Makefile.rd2html ; cd ../


clean.rdoc:
	$(MAKE) clean.all -f Makefile.rdoc


clean.latex2html:
	@for dir in $(SUBDIR) ; do \
	  cd $$dir/ ; \
	  $(MAKE) clean ; \
	  cd ../.. ; \
	done 
