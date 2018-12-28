#= Totorial of deepconv/arare: rakuraku deepconv
#
# Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
# Version::   $Id: namelist.rd,v 1.4 2014/03/04 04:44:05 sugiyama Exp $
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

== その他の設定

以下では, 初期値 (基本場・擾乱場) と物理過程の選択以外の設定項目を一覧する. 

=== 入出力ファイルに関する設定

ファイルに書き込む情報. 

  &fileset_nml
    !ファイルを最終的に変更した人, 組織. 
    FileInstitution = 'GFD Dennou Club (http://www.gfd-dennou.org)'
    ! 実験名 
    FileTitle       = 'cloud moist convection experiment' (デフォルト値)
    ! データファイル作成プログラム名
    FileSource      =  'deepconv/arare5 (http://www.gfd-dennou.org/library/deepconv)' (デフォルト値)
   /

デフォルト値のままで良い場合は変数を書く必要はない. 例えば, 実験名や
プログラム名はデフォルトのままで良い場合には, FileInstitution のみ
指定すれば良い. 

  &fileset_nml
    !ファイルを最終的に変更した人, 組織. 
    FileInstitution = 'XXXX@yyy.org'
   /

=== リスタートファイルに関する設定

リスタートのための入力・出力ファイルを指定する. 

  &restartfileio_nml
    InputFile  = "arare_restart.nc",   ! 入力ファイル名
    OutputFile = "arare_restart2.nc",  ! 出力ファイル名
  /

=== 積分時間に関する設定

時間刻や積分時間, リスタート時間を設定する. 
リスタートする場合の書き方は, ((<こちら|URL:./settings5.htm>)) を参照下さい. 

  &timeset_nml
    DelTimeLong   = 4.0d0      !長いタイムステップ [sec]
    DelTimeShort  = 4.0d-1     !短いタイムステップ(音波関連項) [sec]
    RestartTime   = 0.0d0      !開始時刻 [sec]
    IntegPeriod   = 4800.0d0   !積分時間 [sec]
    DelTimeOutput = 240.0d0    !リスタートファイルに書き出す時間刻み [sec]
  /

=== 格子点数に関する設定

格子点に関する情報を指定します. 
2 次元計算を行う場合の書き方は, ((<こちら|URL:./calc_2d.htm>)) を参照下さい. 
MPI を用いた並列計算のやり方は, ((<こちら|URL:./calc_mpi.htm>)) を参照下さい. 

  &gridset_nml
    xsub  = 1                 ! X 方向の並列数 (デフォルト)
    ysub  = 1                 ! Y 方向の並列数 (デフォルト)
    xdim  = 10                ! X 方向刻み点数 (デフォルト)
    ydim  = 10                ! Y 方向刻み点数 (デフォルト)
    zdim  = 10                ! Z 方向刻み点数 (デフォルト)
    NCMAX = 1                 ! 凝結成分の数 (デフォルト)
    Xmg   = 2                 ! X 方向のマージン (デフォルト)
    Ymg   = 2                 ! Y 方向のマージン (デフォルト)
    Zmg   = 2                 ! Z 方向のマージン (デフォルト)
  /

=== 座標軸に関する設定

座標軸に関する情報を指定します. 

  &axesset_nml
    Xmin  = 0.0d0              ! X 座標の始点 (デフォルト)
    Xmax  = 1.0d4              ! X 座標の終点 (デフォルト)
    Ymin  = 0.0d0              ! X 座標の始点 (デフォルト)
    Ymax  = 1.0d4              ! X 座標の終点 (デフォルト)
    Zmin  = 0.0d0              ! Z 座標の始点 (デフォルト)
    Zmax  = 1.0d4              ! Z 座標の終点 (デフォルト)
  /

デフォルト値は指定しなくても良いです. 

  &axesset_nml
    Xmax  = 1.0d4              ! X 座標の終点 (デフォルト)
    Ymax  = 1.0d4              ! X 座標の終点 (デフォルト)
    Zmax  = 1.0d4              ! Z 座標の終点 (デフォルト)
  /

=== 実験パラメタ(物理・化学)の設定

定圧比熱・定積比熱・分子量・気体定数を指定する方法は以下の 3 通りあります. 
詳しい設定方法については ((<こちら|URL:./settings3.htm>)) を参照下さい. 

* 定圧比熱と分子量を陽に与える場合
* 定圧比熱と気体定数を与える場合
* 物質名とそのモル比を与える場合

デフォルト値は以下の通りです. 

  &constants_nml
    Grav = 9.8d0          !重力 [m/s^2]
    PressBasis = 965.0d0  !温位の基準圧力 [Pa]
    TempSfc = 0.0d0       !地表面温度 [K]
    PressSfc = 0.0d0      !地表面圧力 [Pa]
    TempTop = 0.0d0       !上部境界の温度 [K]
    PressTop = 0.0d0      !上部境界での圧力 [Pa]
    CpDry  = 0.0d0        !乾燥成分の定圧比熱 [J/K kg]
    CpDryMol = 0.0d0      !乾燥成分の定圧比熱 [J/K kg]
    CvDry = 0.0d0         !乾燥成分の定積比熱 [J/K kg]
    MolWtDry = 0.0d0      !乾燥成分の分子量   [kg/mol]
    GasRDry  = 0.0d0      !乾燥成分の気体定数 [J/K kg]
    DayTime = 86400.0d0   ! 1 日の長さ [s]
  /


=== 組成に関する設定

必要に応じて適宜増減させる. 以下は地球大気の場合である.

  &composition_nml
    SpcWetSymbol(1)  = 'H2O-g',  !湿潤成分
    SpcWetSymbol(2)  = 'H2O-s-Cloud', !湿潤成分
    SpcWetSymbol(3)  = 'H2O-s-Rain',  !湿潤成分
    SpcWetMolFr(1)   = 1.0d-2,  !湿潤成分の存在度(モル比)
    SpcWetMolFr(2)   = 0.0d0,  !湿潤成分の存在度(モル比)
    SpcWetMolFr(3)   = 0.0d0,   !湿潤成分の存在度(モル比)
  /


=== 減衰係数の設定(スポンジ層) 

スポンジ層の減衰係数や厚さを設定する. DepthVb は金星実験のように, 
下部境界にスポンジ層を設定する場合に利用する. 

  &damping_nml
    EFTime = 3.0d2             !スポンジ層の減衰係数の e-folding time
    DepthH = 0.0d0             !スポンジ層の厚さ(水平方向)
    DepthV = 0.0d0             !スポンジ層の厚さ(鉛直方向) [上部境界]
    DepthVb= 0.0d0             !スポンジ層の厚さ(鉛直方向) [下部境界] 
  /

=== 減衰係数の設定(音波減衰項)

通常は変更する必要はない.

  &dynamics_nml
   AlphaSound      = 2.0e-7    !音波減衰項の係数
   AlphaNDiff      = 1.0d-4    !数値拡散項の係数
   NDiffRatio      = 1.0d0     !速度に対する粘性を上げる場合は数字を 1 以上にする. 
   beta            = 1.0       !鉛直方向を完全陰解法にする場合には 1.0
                               !クランクニコルソン法にする場合は 0.5 とする. 
   FactorBuoyTemp  = 1.0d0     !浮力 (温度の寄与) の有無
                               !考慮しない場合は値をゼロにする.
   FactorBuoyMolWt = 1.0d0     !浮力 (分子量効果) の有無
                               !考慮しない場合は値をゼロにする.
   FactorBuoyLoading = 1.0d0   !浮力 (荷重効果) の有無
                               !考慮しない場合は値をゼロにする.
  /



=== 出力の全体設定

詳細は((<こちら|URL:settings4.htm>)) をご覧下さい.

  &gtool_historyauto_nml
    FilePrefix = "thermal-moist_" ! 出力ファイルの接頭詞
    IntValue = 10.0,              ! 出力間隔の数値
    IntUnit = 'sec',              ! 出力間隔の単位
  /

=== データ出力の個別設定

出力変数毎の個別の設定を行う. 詳細は((<こちら|URL:settings4.htm>)) をご覧下さい.

  &gtool_historyauto_nml
    Name = ''                     ! 出力変数名
    IntValue = 10.0,              ! 出力間隔の数値
    IntUnit = 'sec',              ! 出力間隔の単位
  !  TimeAverage = .true.,        ! 時間平均の有無
  !  SpaceAverage = .true.,       ! 水平平均の有無
  /


=end
