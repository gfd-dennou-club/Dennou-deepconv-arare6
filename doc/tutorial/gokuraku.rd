#= Totorial of deepconv/arare: rakuraku deepconv
#
# Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
# Version::   $Id: gokuraku.rd,v 1.6 2014/02/26 07:13:10 sugiyama Exp $
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

= ごくらく deepconv/arare5
#= deepconv/arare5 利用の手引
#* 杉山 耕一朗, 小高 正嗣, 山下 達也
#  * 2014/02/05  (杉山 耕一朗) 更新
#  * 2012/05/09  (小高 正嗣) 更新
#  * 2011/03/03  (山下 達也) 更新
#  * 2011/02/24  (山下 達也) 更新
#  * 2009/03/06  (山下 達也) 更新
#  * 2008/06/18  (小高 正嗣) 更新
#  * 2006/11/20  (小高 正嗣) 更新
#  * 2006/10/18  (小高 正嗣) 更新
#  * 2006/09/12  (小高 正嗣) 新規作成

== はじめに

この文書は deeoconv/arare を用いて手軽に実験を行うためのチュートリアル
です. 

#実行プログラムはソースツリーディレクトリ (ここでは /work/deepconv 
#に展開されているとします) 直下の bin ディレクトリ以下にあるとします. 

== deepconv/arare5 のビルド

((<deepconv インストールの手引|URL:../../INSTALL.htm>)) を参考に, 
deepconv/arare5 のビルドを行ってください. 
「インストールの手順」の「ソースのコンパイル」まで行ってください. 
((<いくつかのコンパイラに関する注意書き|URL:./compiler_note.htm>)) も参照下さい. 

ビルドが完了すると, "src/main" ディレクトリ以下に, arare, arare_initdata といった実行ファイルが作成されます. 
また, いくつかのサンプル NAMELIST ファイル (拡張子が .conf のファイル) が "exp_setup_files" に用意されています.


== 実験の実行

下に, いくつかの実験の実行方法についての説明を記します. 始めて deepconv を使う人は, まずどれか 1 つの実験 (例えば「XXXXXXXXXXXXXX」) を実行下さい. 

ちなみに, どの実験も実行するには以下の 4 つのステップで行います.

* 実験ディレクトリの準備
* 初期値の準備
* 実験用データの準備
* 実験の実行

なお, 実験ディレクトリの準備は, 必ずしも必要ありません. ここでは, ある実験の設定や結果が他のものと混ざってしまうことを防ぐために, 各実験ごとにディレクトリを作成しています.

===重力流の実験 (Straka et al., 1993)

重力流の実験を実行する方法を((<こちら|URL:./exp-s93.htm>))で説明します. 

===ケルビンヘルムホルツ不安定の実験 

ケルビンヘルムホルツ不安定の実験を実行する方法を((<こちら|URL:./catseye.htm>))で説明します. 


#* サーマルの上昇 (凝結なし, 乱流: Klemp and Wilhelmson)
#
#   arare-thermal-dry_init-data.conf
#   arare-thermal-dry.conf

#* サーマルの上昇 (凝結: Kessler, H2O のみ, 乱流: Klemp and Wilhelmson)
# 
#   arare-thermal-moist_init-data.conf
#   arare-thermal-moist.conf

#* サーマルの上昇 [木星版] (凝結: Kessler, H2O, NH3, NH4SH, 乱流: Klemp and Wilhelmson)
#
#   arare-jupiter_init-data.conf
#   arare-jupiter.conf

#* サーマルの上昇 [火星版] (凝結: 拡散成長, 乱流: Klemp and Wilhelmson)
#
#   arare-mars_init-data.conf
#   arare-mars.conf

#* Takemi (2007) の初期値
#
#   arare-takemi_init-data.conf

#* Yamasaki (1983) の初期値
#
#   arare-yamasaki_init-data.conf

#* 前バージョン (deepconv/arare4) のヒストリーファイルから初期値生成
#
#   arare-from-arare4_init-data.conf 

#以下, いくつかの実験の実行方法について説明します. 
#
#=== 木星サーマル上昇実験
#
#木星大気を模した条件下で, サーマル上昇実験を実行する方法を
#((<こちら|URL:./jupiter-plumetest.htm>)) で説明します. 
#
#=== 地球サーマル上昇実験
#
#地球大気を模した条件下で, サーマル上昇実験を実行する方法を
#((<こちら|URL:./earth-plumetest.htm>)) で説明します. 
#
#=== 火星サーマル上昇実験
#
#小高他 (2006) で行われた, 火星サーマル上昇実験を実行する方法を
#((<こちら|URL:../../doc-mmc/tutorial/mars-plumetest.htm>)) で説明します. 


== 簡単な解析・可視化

簡単な解析・可視化については, ((<簡単解析・可視化|URL:./visualization.htm>)) を参考にしてください. 

#== 実験条件の変更
##
#実験条件を変更するには namelist を編集する必要があります. 
#namelist ファイルの各項目については 
#((<こちら|URL:namelist.htm>)) で説明します. 

== 実行プログラムの更新

src ディレクトリ以下の編集を行った場合, 
以下のようにすることで実験用ディレクトリの内容を更新することが可能です. 

まず, ソースツリー直下で

 $ make clean

を実行します. 

次に
((<インストールの手引|URL:../../INSTALL.htm>)) の
「インストールの手順」の「ソースのコンパイル」を参考に, 
行ないたい計算に応じてソースのコンパイルを再度行ないます. 

#== 参考文献
#
#* M. Odaka, T. Kitamori, K. Sugiyama, K. Nakajima and Y.-Y. Hayashi, 2006:
#  Numerical simulation of Martian atmospheric convection,
#  Proc. of the 20th ISAS Atmospheric Science Symposium, JAXA/ISAS,
#  103 -- 106


=end
