=begin JA

= 簡単な解析・可視化

# * 森川 靖大 (morikawa), 納多 哲史 (noda), 高橋 芳幸 (yot), 竹広 真一 (takepiro)
# * 杉山 耕一朗 (sugiyama) [copy from dcpam]
#   * $Id: visualization.rd,v 1.2 2014/03/04 14:02:28 sugiyama Exp $

=end JA
=begin EN

= First step analysis and visualization

# * Yasuhiro MORIKAWA (morikawa), Satoshi NODA (noda), Yoshiyuki * O. Takahashi (yot) Shin-ichi Takehiro (takepiro)
# * Ko-ichiro SUGIYAMA (sugiyama) [copy from dcpam]
#   * $Id: visualization.rd,v 1.2 2014/03/04 14:02:28 sugiyama Exp $

=end EN


=begin JA
== 解析・可視化ツールの準備

deepconv/arare5 では入出力するファイルとして
((<Gtool4 NetCDF 規約|URL:http://www.gfd-dennou.org/library/gtool>))
に基づいた NetCDF データを扱います. 

数値実験の結果を解析・可視化するためには, NetCDF
データを取り扱うことのできる解析・可視化ツールが必要です. 
ここでは, ((<電脳 Ruby プロジェクト|URL:http://ruby.gfd-dennou.org/index-j.htm>))
から提供される ((<Gphys|URL:http://www.gfd-dennou.org/library/ruby/products/gphys/>))
を使った可視化の例を紹介します. 

=== 可視化ツールのインストール

((<電脳Ruby謹製品 インストールガイド|URL:http://www.gfd-dennou.org/arch/ruby/tutorial/install/index-j.html>))
を参照してください. 

=end JA

=begin EN
== Preparation of analysis and visualization tools

"deepconv/arare5" input/output NetCDF data based on 
((<Gtool4 NetCDF Conventions|URL:http://www.gfd-dennou.org/library/gtool/index.htm.en>))

Analysis and visualization tools for NetCDF data are
need in order to analyze and visualize results of numerical experiments.
Here, ((<Gphys|URL:http://www.gfd-dennou.org/library/ruby/products/gphys/>))
provided from 
((<Dennou Ruby Project|URL:http://ruby.gfd-dennou.org/index.htm>))
is used via irb. 
Details of usage can be found ((<here|URL:http://www.gfd-dennou.org/library/ruby/tutorial/index-e.html>)).

=== Installation of tool

See ((<Dennou Ruby Products Installation Guide|URL:http://www.gfd-dennou.org/arch/ruby/tutorial/install/>)). 

=end EN

=begin JA

== GPhys/GGraph による解析と可視化

ここでは, 
((<重力流の実験|URL:./exp-s93.htm>))
で得られたデータを GPhys/GGraph を用いて可視化してみることにします. 

まず irb を起動してください. 

  $ irb

以下のような irb のプロンプトが表示されます. 

  irb(main):001:0>

このプロンプトに, 以下のようにコマンドを打ちます. 
左端の数字は行番号で, 打つ必要はありません.

  1: require "numru/ggraph"
  2: include NumRu
  3: gphys = GPhys::IO.open('denscurrent-dry_PTemp.nc', 'PTemp').cut('y'=>0.0)
  4: DCL.gropn(4)
  5: DCL.sgpset('lcntl', false) ; DCL.uzfact(0.7)
  6: GGraph.tone gphys

irb のプロンプトにおいて quit と打つと irb を終了することができます. 

ここでは Temp.nc というファイルの中の Temp という変数を読み込み, 
図示を行っています. 

((<"IMG:exp-s93_img02.png">))

PTemp は x, y, z, t (時間) の 4 次元データですが, 
この実験は水平鉛直 2 次元で行っているので, 
3 行目で y = 0 としています. 
何も指定しないと最後の次元に関しては自動的に 1 番目の要素が選択されます. 
したがってこの図は t=0 での水平鉛直 2 次元での温位を示していることになります

また, 続けて, 下のように時刻を指定することで, 
異なる時刻での最下層の温度分布を描くことができます. 

  7: GGraph.tone gphys.cut('t'=>900)

((<"IMG:exp-s93_img03.png">))

終了させるには 

  8: DCL.grcls
  9: quit

としましょう. 

絵を描くだけでなく, 解析を行うこともできます. 
例として温位とエクスナー関数から温度を計算し, 図示してみましょう. 

  1: require "numru/ggraph"
  2: include NumRu
  3: gphys1 = GPhys::IO.open('denscurrent-dry_PTempAll.nc', 'PTempAll').cut('y'=>0.0)
  4: gphys2 = GPhys::IO.open('denscurrent-dry_ExnerAll.nc', 'ExnerAll').cut('y'=>0.0)
  5: temp = gphys1 * gphys2
  6: DCL.gropn(1)
  7: DCL.sgpset('lcntl', false) ; DCL.uzfact(0.7)
  8: GGraph.tone temp.cut('t'=>900,'x'=>25e3..45e3)
  9: GGraph.contour temp.cut('t'=>900,'x'=>25e3..45e3),false
  10: GGraph.color_bar

((<"IMG:exp-s93_img04.png">))

第 5 行目で温度の計算を温位とエクスナー関数から行っています. 
図示させているのは x = 25 ~ 45 km です. 

このように, GPhys/GGraph を用いると多彩な解析と可視化を実現できます. 
処理が長くなってきたら irb でインタラクティブに行うかわりに, 
エディタを用いてスクリプトファイルを書く方が効率的であり再利用も容易になります. 
より高度な解析・可視化を行う際には, 
((<GPhys チュートリアル|URL:http://ruby.gfd-dennou.org/products/gphys/tutorial/>))
を参照してください. 

=end JA

=begin EN

== Analysis and visualization with GPhys/GGraph

Under construction

=end EN

=begin JA

== GPhys/gpコマンドによる可視化

ここでは, 
((<Polvani et al. (2004) の傾圧不安定波動実験|URL:./exp-p04.htm>))
で得られたデータを GPhys 付属の gp コマンドを用いて可視化してみることにします. 
温度のデータを読み取り図示するには, 

  % gpview denscurrent-dry_PTemp.nc@PTemp,y=0

と入力します. これは denscurrent-dry_PTemp.nc というファイルの中の 
PTemp という変数を読み込み, 図示せよというコマンドです. 

((<"IMG:exp-s93_img05.png">))

PTemp は x, y, z, t (時間) の 4 次元データですが, 
この実験は水平鉛直 2 次元で行っているので, y = 0 としています. 
何も指定しないと最後の次元に関しては自動的に 1 番目の要素が選択されます. 
したがってこの図は t=0 での水平鉛直 2 次元での温位を示していることになります

断面を変えて, 時刻 t=900, x = 25 ~ 50 km にしたければ, 

  % gpview denscurrent-dry_PTemp.nc@PTemp,y=0,x=25e3:50e3,t=900 --nocont 

と, カンマで区切って断面を指定することができます. 

((<"IMG:exp-s93_img06.png">))

アニメーションも簡単に見ることができます. 

  % gpview denscurrent-dry_PTemp.nc@PTemp,y=0,x=25e3:50e3 --anim t --nocont

座標を 3 つ指定して 1 次元データにすると折れ線グラフを描けます. 

  % gpview denscurrent-dry_PTemp.nc@PTemp,y=0,x=25e3:50e3,t=900 --exch

--exch オプションは縦軸と横軸を入れ換える操作を指示するものです. 

((<"IMG:exp-s93_img07.png">))

gpview には他にもいろいろなオプションがあります. 
gpview --help とするとオプションと使い方の例が表示されます. 
gp コマンドシリーズには他にもいろいろなものが用意されています. 
主なものは以下の通りです. 

: gpvect
  2 次元ベクトル図の表示

: gpprint
  データの数値の出力表示

: gplist
  ファイルに格納されている変数のリストを表示

: gpmaxmin
  データの最大・最小値の表示

gp コマンドシリーズは 1 行入力ですぐに結果を表示できるのが特徴です. 
そのためクイックルックや計算の途中でのデータチェック等に便利です. 
しかしながら複数の変数を組み合わせた解析や可視化はできません. 
本格的なデータ解析と可視化を行うには, 
((<"GPhys/GGraph による解析と可視化">))が適しているでしょう. 

=end JA

=begin EN

== Visualization with GPhys/gpcommands

Under construction

=end EN

=begin HTML
<hr />
<small>
  $Id: visualization.rd,v 1.2 2014/03/04 14:02:28 sugiyama Exp $
</small>
=end HTML

