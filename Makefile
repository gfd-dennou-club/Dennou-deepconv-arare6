#= Makefile for deepconv/arare source tree.
#
# Authors::   Masatsugu ODAKA
# Version::   $Id: Makefile,v 1.3 2014/02/21 04:19:37 deepconv Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:COPYRIGHT]
#
#
#== History
#
# * 2007/10/19 (ODAKA Masatsugu): Update
# * 2006/09/28 (ODAKA Masatsugu): Update
# * 2006/09/20 (ODAKA Masatsugu): Update
# * 2006/09/13 (ODAKA Masatsugu): Update
# * 2006/09/12 (ODAKA Masatsugu): Update
# * 2006/03/19 (SUGIYAMA Ko-ichiro):
# * 2005/04/26 (ODAKA Masatsugu): 
# * 2005/04/22 (SUGIYAMA Ko-ichiro):
# * 2005/01/27 (ODAKA Masatsugu): Update
# * 2003/02/05 (SUGIYAMA Ko-ichiro):
#
########################################################################
#== Settings
#
# [JAPANESE] 各種設定項目

INCLUDEFILE = Config.mk
include $(INCLUDEFILE)  # Include file
                        # インクルードファイル

GuideFiles = INSTALL CREDITS
                      # Documentation located in top directory       (optional)
                      # トップディレクトリに置くドキュメントファイル (任意)

# End Settings  [JA] 設定項目ここまで
########################################################################

######################################################################
#== Set GuideRDFiles
GuideRDFiles = $(GuideFiles:=.rd)
GuideJA      = $(GuideFiles:=.htm)
GuideEN      = $(GuideFiles:=.htm.en)

#== Set subdirectories
#DOCSUBDIR = $(SRCDIR) $(DOCDIR) 

######################################################################
#== Rules
#


all: 
	cd $(SRCDIR) ; $(MAKE) all

mpi: 
	cd $(SRCDIR) ; $(MAKE) mpi

3d:
	cd $(SRCDIR) ; $(MAKE) 3d

install:
	install -d $(DEST_INC)
	install -d $(DEST_BIN)
	install -d $(DEST_LIB)
	install $(MODDIR)/*.mod $(DEST_INC)
	install $(LIBDIR)/lib$(LIBNAME).a $(DEST_LIB)
	install $(BINDIR)/* $(DEST_BIN)

doc: rd2html
	@for dir in $(DOCDIR) ; do \
	  cd $$dir ; \
	  $(MAKE) doc ; \
	  cd ../ ; \
	done


rd2html: $(GuideFiles)
	@for file in $^ ; do \
	  echo $(CP) $${file} $${file}.rd ;\
	  $(CP) $${file} $${file}.rd ;\
	done
	rake
#	$(MAKE) -f Makefile.rd2html
#	$(RM) $(GuideRDFiles)


latex2html: 
	@for dir in $(DOCDIR) ; do \
	  cd $$dir ; \
	  $(MAKE) latex2html ; \
	  cd ../ ; \
	done


clean.all: clean.doc clean.latex2html clean.config clean
	-$(RM) -f Config.mk chkfort.cfg  config.log  config.status 

clean.doc: clean.rd2html
	cd $(DOCDIR) ; $(MAKE) clean ; cd ../

clean.rd2html:
	-$(RM) -f $(GuideJA) $(GuideEN)


clean.latex2html:
	cd $(DOCDIR) ; $(MAKE) clean.latex2html ; cd ../


clean.config:
	-$(RM) -f chkfort.cfg config.cache config.log config.status

clean:
	cd $(SRCDIR); $(MAKE) clean ; cd ../
	-$(RM) -r $(LIBDIR)
	-$(RM) -r $(MODDIR)
	-$(RM) -r $(BINDIR)
