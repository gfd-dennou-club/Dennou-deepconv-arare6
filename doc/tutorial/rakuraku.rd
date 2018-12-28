#= Totorial of deepconv/arare: rakuraku deepconv
#
# Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
# Version::   $Id: rakuraku.rd,v 1.5 2014/03/04 04:44:05 sugiyama Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:../../COPYRIGHT]
#
#
########################################################################
#
=begin TOPLINK
[((<GFD Dennou Club|URL:http://www.gfd-dennou.org>)) |
((<deepconv Project|URL:http://www.gfd-dennou.org/library/deepconv>))]
[((<deepconv Reference Manual|URL:../index.htm>))]
=end TOPLINK

=begin

= らくらく deepconv/arare5
#* 杉山 耕一朗, 小高 正嗣, 山下 達也
#  * 2012/05/09  (小高 正嗣) 更新
#  * 2009/03/06  (山下 達也) 更新
#  * 2008/06/18  (小高 正嗣) 更新
#  * 2006/11/20  (小高 正嗣) 更新
#  * 2006/10/18  (小高 正嗣) 更新
#  * 2006/09/12  (小高 正嗣) 新規作成


== はじめに

この文書は deeoconv/arare の利用の手引である. ここでは ((<インストール
の手引|URL:../../INSTALL.htm>)) により, すでにソースコードのコンパイル
と実行プログラムのインストールは終了しているものとし, 自分で基本場と初
期条件を設定し実行する方法について簡単に解説する. デフォルトの設定を用
いて実行する方法については((<ごくらく deepconv|URL:gokuraku.htm>))を参
照されたい.

== Contents 

*設定ファイルを用いた実験設定の変更
  * ((<解像度を変更するには|URL:./settings1.htm>))
  * ((<時間刻み・積分期間を変更するには|URL:./settings2.htm>))
  * ((<物理定数・惑星に関する定数を変更するには|URL:./settings3.htm>))
  * ((<出力設定を変更するには|URL:./settings4.htm>))
  * ((<リスタート計算を行うには|URL:./settings5.htm>))

#* ソースの変更・追加

* ((<並列計算を行うには|URL:./calc_mpi.htm>))
  * deepconv/arare5 は MPI による並列化のみ行っている. OpenMP は用いていない. 

* ((<2 次元の計算を行うには|URL:./calc_2d.htm>))

* ((<凝結物の扱いに関する解説|URL:./condense.htm>))

#* ((<初期値生成用設定ファイル (NAMELIST)|URL:./namelist_init.htm>))

* 設定ファイル (NAMELIST)
  * ((<初期値 (基本場・擾乱場) の選択|URL:./namelist_init.htm>))
  * ((<物理過程の選択|URL:./namelist_main.htm>))
  * ((<その他の設定|URL:./namelist.htm>))

#* 付属スクリプトの解説

* ((<変数名・微分平均演算について|URL:grid_operation.htm>))

* テスト計算について [準備中]
  * 飽和蒸気圧のチェック, 微分演算の精度チェック, など. 


=end 
