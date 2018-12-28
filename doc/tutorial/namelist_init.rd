#= Totorial of deepconv/arare: rakuraku deepconv
#
# Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
# Version::   $Id: namelist_init.rd,v 1.1 2014/03/04 04:44:05 sugiyama Exp $
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

= 設定ファイルの概要


== 初期場の設定

基本場 (FlagBasic) と擾乱場 (FlagDisturb...) を設定する. 
基本場の設定は必須. 擾乱場は 1 つ以上設定する必要がある. 

  &initialdata_nml
    !
    ! [必須] 基本場の設定. 以下から選択すること. 
    !
    ! FlagBasic = "IsoThermal"   ! 計算領域の温度一定
    ! FlagBasic = "Dry"          ! 乾燥断熱な基本場
    ! FlagBasic = "Sounding"     ! サウンディングファイルより与える
    ! FlagBasic = "Yamasaki1983" ! Yamasaki et al (1983) を模した基本場
    ! FlagBasic = "Takemi2007"   ! Takemi (2007) を模した基本場
    ! FlagBasic = "Toon2002"     ! Toon et al. (2002) を模した基本場
    ! FlagBasic = "Baker1998"    ! deepconv/arare4 の出力を利用
    ! FlagBasic = "Arare4"       ! deepconv/arare4 の出力を利用
    ! FlagBasic = "Arare4mmc"    ! deepconv/arare4 (火星計算) の出力を利用
    ! 
    ! 温位の擾乱場. 以下を選択可能. 
    !
    ! FlagDisturbPTemp = "GaussXY"   ! ガウシアン (XY 平面, Z 方向一定)
    ! FlagDisturbPTemp = "GaussXZ"   ! ガウシアン (XZ 平面, Y 方向一定)
    ! FlagDisturbPTemp = "GaussYZ"   ! ガウシアン (YZ 平面, X 方向一定)
    ! FlagDisturbPTemp = "GaussXYZ"  ! ガウシアン (XY 平面, 3 次元的な分布)
    ! FlagDisturbPTemp = "Random"    ! ランダム
    ! FlagDisturbPTemp = "Rectangle" ! 矩形
    ! FlagDisturbPTemp = "CosXY"     ! cos (XY 平面, Z 方向一定)
    ! FlagDisturbPTemp = "CosXZ"     ! cos (XZ 平面, Y 方向一定)
    ! FlagDisturbPTemp = "CosYZ"     ! cos (YZ 平面, X 方向一定)
    ! FlagDisturbPTemp = "CosXYZ"    ! cos 
    ! FlagDisturbPTemp = "ConeXY"    ! 円錐
    ! FlagDisturbPTemp = "ConeXZ"    ! 円錐
    ! FlagDisturbPTemp = "ConeYZ"    ! 円錐
    ! FlagDisturbPTemp = "tanh"      ! tanh 型のシア
    !
    ! エクスナー関数の擾乱場. 以下を選択可能. 
    ! 
    ! FlagDisturbExner = "GaussXY"   ! ガウシアン (XY 平面, Z 方向一定)
    ! FlagDisturbExner = "GaussXZ"   ! ガウシアン (XZ 平面, Y 方向一定)
    ! FlagDisturbExner = "GaussYZ"   ! ガウシアン (YZ 平面, X 方向一定)
    ! FlagDisturbExner = "GaussXYZ"  ! ガウシアン 
    ! FlagDisturbExner = "CosXY"     ! cos (XY 平面, Z 方向一定)
    ! FlagDisturbExner = "CosXZ"     ! cos (XZ 平面, Y 方向一定)
    ! FlagDisturbExner = "CosYZ"     ! cos (YZ 平面, X 方向一定)
    ! FlagDisturbExner = "CosXYZ"    ! cos
    ! FlagDisturbExner = "ConeXY"    ! 円錐
    ! FlagDisturbExner = "ConeXZ"    ! 円錐
    ! FlagDisturbExner = "ConeYZ"    ! 円錐
    !
    ! 凝結成分の混合比. 以下を選択可能. 
    !
    ! FlagDisturbQMix = "GaussXY"    ! ガウシアン (XY 平面, Z 方向一定)
    ! FlagDisturbQMix = "GaussXZ"    ! ガウシアン (XZ 平面, Y 方向一定)
    ! FlagDisturbQMix = "GaussYZ"    ! ガウシアン (YZ 平面, X 方向一定)
    ! FlagDisturbQMix = "GaussXYZ"   ! ガウシアン
    ! FlagDisturbQMix = "Dryreg"     ! 乾燥領域を設定
#    ! FlagDisturbQMix = "Moist"      ! ?
    ! FlagDisturbQMix = "CosXY"      ! cos (XY 平面, Z 方向一定)
    ! FlagDisturbQMix = "CosXZ"      ! cos (XZ 平面, Y 方向一定)
    ! FlagDisturbQMix = "CosYZ"      ! cos (YZ 平面, X 方向一定)
    ! FlagDisturbQMix = "CosXYZ"     ! cos
    ! FlagDisturbQMix = "ConeXY"     ! 円錐
    ! FlagDisturbQMix = "ConeXZ"     ! 円錐
    ! FlagDisturbQMix = "ConeYZ"     ! 円錐
    !
    ! 速度場の混合比. 以下を選択可能. 
    !
    ! FlagDisturbWind = "Takemi2007" ! Takemi (2007) を模した速度場
    ! FlagDisturbWind = "Sounding"   ! サウンディングファイルから速度場を与える
    ! FlagDisturbWind = "Zonal"      ! 水平一様な速度場
    ! FlagDisturbWind = "tanh"       ! tanh 型のシア
    !
    ! FlagDisturb = "Arare4"         ! deepconv/arare4 の出力を利用
    ! FlagDisturb = "Arare4mmc"      ! deepconv/arare4 (火星計算) の出力を利用
  /

== 基本場の設定の際に必要となる情報

=== FlagBasic = "IsoThermal" の場合

特になし. 

=== FlagBasic = "Dry" の場合

  &initialdata_basic_nml
     Humidity
     TempTr
     DHeight
     HeightTr
  /

=== FlagBasic = "Sounding"  の場合

  &initialdata_sounding_nml
     SoundingFile = ""     ! サウンディングファイル
     AltCol       = 0      ! 「高度」の列番号 (サウンディングファイル内)
     TempCol      = 0      ! 「温度」の列番号 (サウンディングファイル内)
     PressCol     = 0      ! 「圧力」の列番号 (サウンディングファイル内)
     VelXCol      = 0      ! 「X 方向の速度」の列番号 (サウンディングファイル内)
     VelYCol      = 0      ! 「Y 方向の速度」の列番号 (サウンディングファイル内)
     !
     ! 以下の情報は, 指定した高度より上空を温度一定の成層圏に設定する場合に必要.
     !
     AltTr        = 0      ! 対流圏界面の高度
     DelAlt       = 4.0d3  ! 温度減率をゼロに緩和するまでの距離を指定
  /

=== FlagBasic = "Yamasaki1983" 

特に設定する必要はない. 
具体的な値は, ((<初期値の計算例|URL:http://www.gfd-dennou.org/library/deepconv/sample/2011-06-23_sugiyama/Init/exp.htm>)) 参照.


=== FlagBasic = "Takemi2007" 

基本場設定する. 
具体的な値は, ((<初期値の計算例|URL:http://www.gfd-dennou.org/library/deepconv/sample/2011-06-23_sugiyama/Init/exp.htm>)) 参照.

  &initialdata_takemi2007_nml
    ! 
    ! 基本場の選択
    !
    ! FlagEnv = "MidLat_Q10"
    ! FlagEnv = "MidLat_Q12"
    ! FlagEnv = "MidLat_Q14"
    ! FlagEnv = "MidLat_Q16"
    ! FlagEnv = "MidLat_Q16DRY1"
    ! FlagEnv = "MidLat_Q16DRY2"
    ! FlagEnv = "MidLat_Q18"
    ! FlagEnv = "Tropic_Q18"
    ! FlagEnv = "Tropic_Q18DRY1"
    ! FlagEnv = "Tropic_Q18DRY2"
    ! FlagEnv = "Tropic_Q18DRY3"
    !
    !! 速度場 (鉛直プロファイル) の選択
    !!
    !! FlagWind = "LowLevel"
    !! FlagWind = "MiddleLevel"
    !! FlagWind = "HighLevel"
    !!
    !! 地表面付近の速度
    !!
    !! VelXSfc = 0.0d0
  /

=== FlagBasic = "Toon2002"

特に設定する必要はない. 
具体的な値は, ((<初期値の計算例|URL:http://www.gfd-dennou.org/library/deepconv/sample/2011-06-23_sugiyama/Init/exp.htm>)) 参照.

=== FlagBasic = "Baker1998"

特に設定する必要はない. 

=== FlagBasic = "Arare4" or "Arare4mmc" の場合

deepconv/arare4 の出力を元に基本場・擾乱場を決める場合には以下の設定を行う. 

  &arare4fileio_nml
    Arare4Prefix = ""   ! deepconv/arare4 のファイルの接頭詞を指定
  /


== 擾乱場の設定に必要となる情報 

=== ガウシアン [温位, エクスナー関数, 混合比]

FlagDisturbPTemp, FlagDisturbExner, FlagDisturbQMix などで, 
(({"Gauss..."})) を選択した場合は, 
initialdata_disturb_gauss_nml の設定を行う. 

  &initialdata_disturb_gauss_nml
     PTempMax = 1.0d0      ! 温位擾乱の最大値
     ExnerMax = 0.01d0     ! エクスナー関数の最大値
     QMixMax  = 1.0d-2     ! 混合比 (気体) の最大値
     Xc       = 1000.0d0   ! 座標の中心 (X 方向)
     Xr       = 100.0d0    ! 半径 (X 方向)
     Yc       = 1000.0d0   ! 座標の中心 (Y 方向)
     Yr       = 100.0d0    ! 半径 (Y 方向)
     Zc       = 1000.0d0   ! 座標の中心 (Z 方向)
     Zr       = 100.0d0    ! 半径 (Z 方向)
  /

注: 必要な値のみ設定すれば良い. 例えば, (({FlagDisturbPTemp = "GaussXY"})) を
指定した場合には, 以下のように設定すれば十分である. 

  &initialdata_disturb_gauss_nml
     PTempMax = 1.0d0      ! 温位擾乱の最大値
     Xc       = 1000.0d0   ! 座標の中心 (X 方向)
     Xr       = 100.0d0    ! 半径 (X 方向)
     Yc       = 1000.0d0   ! 座標の中心 (Y 方向)
     Yr       = 100.0d0    ! 半径 (Y 方向)
  /


=== 三角関数 [温位, エクスナー関数, 混合比]

FlagDisturbPTemp, FlagDisturbExner, FlagDisturbQMix などで, 
(({"Cos..."})) を選択した場合は, 
initialdata_disturb_cos_nml の設定を行う. 

  &initialdata_disturb_cos_nml
     PTempMax = 1.0d0      ! 温位擾乱の最大値
     ExnerMax = 0.01d0     ! エクスナー関数の最大値
     QMixMax  = 1.0d-2     ! 混合比 (気体) の最大値
     Xc       = 1000.0d0   ! 座標の中心 (X 方向)
     Xr       = 100.0d0    ! 半径 (X 方向)
     Yc       = 1000.0d0   ! 座標の中心 (Y 方向)
     Yr       = 100.0d0    ! 半径 (Y 方向)
     Zc       = 1000.0d0   ! 座標の中心 (Z 方向)
     Zr       = 100.0d0    ! 半径 (Z 方向)
  /

注: 必要な値のみ設定すれば良い. 例えば, (({FlagDisturbPTemp = "CosXY"})) を
指定した場合には, 以下のように設定すれば十分である. 

  &initialdata_disturb_cos_nml
     PTempMax = 1.0d0      ! 温位擾乱の最大値
     Xc       = 1000.0d0   ! 座標の中心 (X 方向)
     Xr       = 100.0d0    ! 半径 (X 方向)
     Yc       = 1000.0d0   ! 座標の中心 (Y 方向)
     Yr       = 100.0d0    ! 半径 (Y 方向)
  /


=== 円錐 [温位, エクスナー関数, 混合比]

FlagDisturbPTemp, FlagDisturbExner, FlagDisturbQMix などで, 
(({"Cone..."})) を選択した場合は, 
initialdata_disturb_cone_nml の設定を行う. 

  &initialdata_disturb_cone_nml
     PTempMax = 1.0d0      ! 温位擾乱の最大値
     ExnerMax = 0.01d0     ! エクスナー関数の最大値
     QMixMax  = 1.0d-2     ! 混合比 (気体) の最大値
     Xc       = 1000.0d0   ! 座標の中心 (X 方向)
     Xr       = 100.0d0    ! 半径 (X 方向)
     Yc       = 1000.0d0   ! 座標の中心 (Y 方向)
     Yr       = 100.0d0    ! 半径 (Y 方向)
     Zc       = 1000.0d0   ! 座標の中心 (Z 方向)
     Zr       = 100.0d0    ! 半径 (Z 方向)
  /

注: 必要な値のみ設定すれば良い. 例えば, (({FlagDisturbPTemp = "ConeXY"})) を
指定した場合には, 以下のように設定すれば十分である. 

  &initialdata_disturb_cone_nml
     PTempMax = 1.0d0      ! 温位擾乱の最大値
     Xc       = 1000.0d0   ! 座標の中心 (X 方向)
     Xr       = 100.0d0    ! 半径 (X 方向)
     Yc       = 1000.0d0   ! 座標の中心 (Y 方向)
     Yr       = 100.0d0    ! 半径 (Y 方向)
  /

=== ランダム [温位]

とある高度にランダムな擾乱を設定する場合 (FlagDisturbPTemp = "Random"), 
initialdata_disturb_random_nml の設定を行う. 

  &initialdata_disturb_random_nml
     PTempMax = 1.0d0      ! 擾乱の最大値
     Zpos     = 1.0d3      ! 擾乱を置く高度
  /

=== 矩形 [温位], DryReg [混合比]

矩形領域のみ温度を変化させる場合 (FlagDisturbPTemp = "Rectangle"), 
initialdata_disturb_rectangle_nml の設定を行う. 

  &initialdata_disturb_rectangle_nml
     PTempMax = 1.0d0      ! 温位の増加分
     XposMin  = 0.0d0      ! X 方向の始点
     XposMax  = 1.0d2      ! X 方向の終点
     YposMin  = 0.0d0      ! Y 方向の始点
     YposMax  = 1.0d2      ! Y 方向の終点
     ZposMin  = 0.0d0      ! Z 方向の始点
     ZposMax  = 1.0d2      ! Z 方向の終点
  /

矩形領域の凝結物の混合比をゼロにする場合 (FlagDisturbQMix = "DryReg") も, 
initialdata_disturb_rectangle_nml の設定を行う. 

  &initialdata_disturb_rectangle_nml
     XposMin  = 0.0d0      ! X 方向の始点
     XposMax  = 1.0d2      ! X 方向の終点
     YposMin  = 0.0d0      ! Y 方向の始点
     YposMax  = 1.0d2      ! Y 方向の終点
     ZposMin  = 0.0d0      ! Z 方向の始点
     ZposMax  = 1.0d2      ! Z 方向の終点
  /

=== tanh [温位, 速度場] 

ケルビン・ヘルムホルツ不安定の初期値で用いるような tanh 型のシアを与える場合 
(FlagDisturbPTemp = "tanh" & FlagDisturbWind = "tanh"), 
initialdata_disturb_tanh_nml を設定する. 

  &initialdata_disturb_tanh_nml
     PTempMean  = 平均温位
     PTempDel   = シア層の上下の温位差     
     VelMean    = 平均風速
     Zc         = シア層の中心の高度
     Zr         = シア層の厚さ
  /

=== Zonal [速度場]

一様な風速場を与える場合 (FlagDisturbWind = "Zonal"), 
initialdata_disturb_wind_zonal_nml を設定する. 

  &initialdata_disturb_wind_zonal_nml/ &
    VelX0 = 10.0d0
    VelY0 = 10.0d0
    VelZ0 = 10.0d0
  /

=== Sounding [速度場]

  &initialdata_sounding_nml
     SoundingFile = ""     ! サウンディングファイル
     VelXCol      = 0      ! 「X 方向の速度」の列番号 (サウンディングファイル内)
     VelYCol      = 0      ! 「Y 方向の速度」の列番号 (サウンディングファイル内)
  /

=== Takemi2007 [速度場]

速度場を設定する. 

  &initialdata_takemi2007_nml
    !! 
    !! 基本場の選択
    !!
    !! FlagEnv = "MidLat_Q10"
    !! FlagEnv = "MidLat_Q12"
    !! FlagEnv = "MidLat_Q14"
    !! FlagEnv = "MidLat_Q16"
    !! FlagEnv = "MidLat_Q16DRY1"
    !! FlagEnv = "MidLat_Q16DRY2"
    !! FlagEnv = "MidLat_Q18"
    !! FlagEnv = "Tropic_Q18"
    !! FlagEnv = "Tropic_Q18DRY1"
    !! FlagEnv = "Tropic_Q18DRY2"
    !! FlagEnv = "Tropic_Q18DRY3"
    !!
    ! 速度場 (鉛直プロファイル) の選択
    !
    ! FlagWind = "LowLevel"
    ! FlagWind = "MiddleLevel"
    ! FlagWind = "HighLevel"
    !
    ! 地表面付近の速度
    !
    ! VelXSfc = 0.0d0
  /

=== "Arare4" or "Arare4mmc" の場合

deepconv/arare4 の出力を元に基本場・擾乱場を決める場合には以下の設定を行う. 

  &arare4fileio_nml
    Arare4Prefix = ""   ! deepconv/arare4 のファイルの接頭詞を指定
  /

=end
