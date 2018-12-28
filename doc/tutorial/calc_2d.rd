=begin JA

= 水平鉛直 2 次元の計算を行うには

# * 杉山耕一朗 (sugiyama)
#   * $Id: calc_2d.rd,v 1.1 2014/03/01 20:11:21 sugiyama Exp $


本文書では水平鉛直 2 次元の計算を実行するための方法を示す. 


== 設定ファイルの書き方

設定ファイル(NAMELIST ファイル) の &gridset_nml, &axesset_nml を変更す
るだけである. (({ydim = 1})) とすることで, deepconv/arare5 は水平鉛直
2 次元のモデルとして動作する. 

    &gridset_nml
      xdim  = 50                 ! X 方向刻み点数
      ydim  = 1                  ! Y 方向刻み点数
      zdim  = 50                 ! Z 方向刻み点数
      NCMAX = 1                  ! 化学種の数
    /

    &axesset_nml
      Xmax  = 1.0d3             ! X 座標の終点
      Zmax  = 6.5d2             ! Z 座標の終点
    /

このとき, Ymin, Ymax は指定する必要は無い. 


=end JA


=begin HTML
<hr />
<small>
  $Id: calc_2d.rd,v 1.1 2014/03/01 20:11:21 sugiyama Exp $
</small>
=end HTML

