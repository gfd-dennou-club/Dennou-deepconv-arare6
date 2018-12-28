# 表題	deepconv 離散化マニュアル用 Makefile
# Copyright (C) GFD_Dennou_Club, 2000.  All rights reserved.
#
# 履歴  1999/07/21  小高正嗣, 林祥介, 
#       1999/11/20  竹広真一
#       2000/01/15  杉山耕一朗
#       2000/04/11  竹広真一
#       2011/12/06  小高正嗣: Mkinclude ファイルと統合
#       2014/02/16  竹広真一
#
##################################################################
#
# tex ソースのメインファイル名
#
TARGET=main
#
##################################################################
#
# 現在の作業ディレクトリ
#
WORKDIR=$(PTOPDIR)/src 
#
# 生成物のインストール先
#
TARGETDIR=../pub
TARGETHTMDIR=../htm
##################################################################
#
# コマンド定義
#
PLATEX= platex -kanji=euc
### PLATEX= ajlatex
DVIPS = dvips -f 
### DVIPS = dvi2ps 
PS2PDF= ps2pdf
DVI2PDF= dvipdfmx
MKDIR = mkdir
RM    = rm -f
LATEX2HTML = latex2html
LATEX2HTMLOPTIONS = -no_footnode -local_icons -short_extn \
		-bottom_navigation -no_contents_in_navigation -no_auto_link \
		-up_url $(HREFTOPDIR)/index.htm -up_title $(HREFTOPDIRNAME) \
		-split 3 -show_section_numbers -dir $(TARGETHTMDIR) 

##################################################################
#
# rules and dependencies 
#

#all: installdir pdf latex2html
all: installdir pdf

installdir:
	if [ ! -d $(TARGETDIR) ] ; then \
	$(MKDIR) $(TARGETDIR) ; \
	fi ;\

$(TARGET).dvi: $(TARGET).tex
	$(PLATEX) $(TARGET).tex
	$(PLATEX) $(TARGET).tex

$(TARGET).ps: $(TARGET).dvi
	$(DVIPS) $(TARGET).dvi > $(TARGET).ps

pdf: $(TARGETDIR)/$(TARGET).pdf

$(TARGETDIR)/$(TARGET).pdf: $(TARGET).dvi
	$(DVI2PDF) $(TARGET).dvi 
	mv $(TARGET).pdf $(TARGETDIR)/

gif: $(TARGET).ps
	$(DCNOTE-PS2GIF) -d $(TARGETDIR) $(TARGET).ps

htm: gif
	$(DCNOTE-HTMS) $(TARGETDIR)

latex2html: $(TARGET).dvi
	if [ ! -d $(TARGETHTMDIR) ] ; then \
	$(MKDIR) $(TARGETHTMDIR) ; \
	fi ;\
	$(LATEX2HTML) $(LATEX2HTMLOPTIONS) $(TARGET).tex

clean: 
	$(RM) *~ *.bak *.aux *.log *.toc $(TARGET).ps $(TARGET).dvi

clean.all: clean
	@test -d $(TARGETDIR) && eval "cd $(TARGETDIR) ; $(RM) *.gif *.htm *.pdf" || true
	@test -d $(TARGETDIR) && rmdir $(TARGETDIR) || true
	@test -d $(TARGETHTMDIR) && eval "cd $(TARGETHTMDIR) ; $(RM) *.css *.htm *.log *.pl *.png *.old images.aux images.tex images.idx WARNINGS" || true
	@test -d $(TARGETHTMDIR) && rmdir $(TARGETHTMDIR) || true
