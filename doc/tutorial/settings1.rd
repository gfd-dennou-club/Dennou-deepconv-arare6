=begin JA

= 設定ファイルを用いた実験設定の変更

# * 杉山耕一朗 (sugiyama)
#   * $Id: settings1.rd,v 1.1 2014/03/01 20:11:22 sugiyama Exp $


本文書では設定ファイル(NAMELIST ファイル) を用いた実験設定の変更方法に
ついて記す. 

設定ファイルを変更した後の実際の計算実行の方法については「ごくらくdcpam5」
http://www.gfd-dennou.org/library/dcpam/dcpam5/dcpam5_latest/doc/gokuraku/
を参照されたい.

== 解像度を変更するには

解像度は, 設定ファイル(NAMELIST ファイル) に &gridset_nml,
&axesset_nml を用いて設定する. 

例えば設定ファイルに以下のように書かれている場合, x, y, z 方向の格子サ
イズは 50, 50, 50 であり, dx, dy, dz はそれぞれ 20 m (1000 / 50), 20 m, 
13 m (650 / 50) である.

    &gridset_nml
      xdim  = 50                 ! X 方向刻み点数
      ydim  = 50                 ! Y 方向刻み点数
      zdim  = 50                 ! Z 方向刻み点数
      NCMAX = 1                  ! 化学種の数
    /

    &axesset_nml
      Xmax  = 1.0d3             ! X 座標の終点
      Ymax  = 1.0d3             ! Y 座標の終点
      Zmax  = 6.5d2             ! Z 座標の終点
    /

解像度を 1/2 倍する場合は, 以下のように修正する. この場合, x, y, z 方向
の格子サイズは 100, 1, 100 であり, dx, dy, dz はそれぞれ, 10 m (1000 /
100), 10 m, 6.5 m (650 / 100) である.

    &gridset_nml
      xdim  = 100                ! X 方向刻み点数
      ydim  = 100                ! Y 方向刻み点数
      zdim  = 100                ! Z 方向刻み点数
      NCMAX = 1                  ! 化学種の数
    /

    &axesset_nml
      Xmax  = 1.0d3             ! X 座標の終点
      Ymax  = 1.0d3             ! Y 座標の終点
      Zmax  = 6.5d2             ! Z 座標の終点
    /

なお, 座標軸の下限は何も指定しない場合はゼロである. 例えば, 高度 350 m から
高度 650 m を計算する場合は, 

    &axesset_nml
      Xmax  = 1.0d3             ! X 座標の終点
      Ymax  = 1.0d3             ! Y 座標の終点
      Zmin  = 3.5d2             ! Z 座標の始点
      Zmax  = 6.5d2             ! Z 座標の終点
    /

とする. 




=end JA


=begin HTML
<hr />
<small>
  $Id: settings1.rd,v 1.1 2014/03/01 20:11:22 sugiyama Exp $
</small>
=end HTML

