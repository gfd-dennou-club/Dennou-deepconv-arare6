#= Totorial of deepconv/arare: rakuraku deepconv
#
# Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
# Version::   $Id: namelist_main.rd,v 1.1 2014/03/04 04:44:05 sugiyama Exp $
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


== 物理過程の選択

deepconv_main_nml において, 数値実験に利用する物理過程を選択する. 

  &deepconv_main_nml
    !
    ! 乱流過程
    !FlagTurbMethod  = "KW1978",
    !  KW1978 (Klemp and Wilhelmson 1978)
    !  ConstKm (乱流拡散係数一定)
    !
    ! 雲微物理過程
    !FlagCloudMethod = "K1969",
    !   K1969 (暖かい雨のパラメタリゼーション (Kessler, 1969))
    !   MarsCond (拡散成長)
    !
    ! 放射過程
    !FlagRadMethod = "HeatConst"
    !  HeatConst (指定した高度を一様冷却)
    !  HeatVary  (指定した高度を一様冷却, 冷却率は高度の関数)
    !  HeatBalance (系に与える加熱冷却が釣り合うように熱強制を与える)
    !  Baker1998 (金星用, Baker et al., 1998)
    !  Sounding  (ファイルから加熱率・冷却率を与える)
    !
    ! 地表面過程
    !FlagSurfaceMethod = "Diff"
    !  Diff (境界から拡散的に熱と物質を与える, 下部境界の温度・物質量固定)
    !  Bulk (バルク法)
    !  Const (上部, 下部境界からの一定の熱・運動量・物質フラックスを与える)
    !  Baker1998 (金星用, Baker et al., 1998)
    !
    ! デバッグ用
    !FlagDebugMethod = "Const"
    !  WindConst (速度一定, 与えられた速度で計算結果を上書きする.)
    !  NoTendencyLong (長い時間ステップで計算された tendency をゼロとする.)
    !    
   /    


== 乱流過程

=== FlagTurbMethod  = "KW1978",

Klemp and Wilhelmson (1978) の 1.5 次のクロージャーに基づいて乱流過程を計算する場合には
以下の値を設定することが出来る. デフォルト値のままで良い場合は, 値を設定しなくて良い. 

  &turbulence_kw1978_nml
     Cm     = 2.0d-1          !乱流エネルギー診断式の係数
     KmMax  = 800.0d0         !乱流拡散係数の最大値
     FlagDExnerDtTurb =.true. !圧力方程式に乱流拡散項を考慮するかのスイッチ
  /

=== FlagTurbMethod  = "ConstKm"

乱流拡散係数の値を一定値にする場合には, 以下の値を設定することが出来る. 
デフォルト値のままで良い場合は, 値を設定しなくて良い. 

  &turbulence_constKm_nml
     Cm     = 2.0d-1          !乱流エネルギー診断式の係数
     MixLen = 0.0d0           !平均混合距離
     ConstKm = 0.0d0          !運動量に対する乱流拡散係数
     ConstKh = 0.0d0          !熱に対する乱流拡散係数
     FlagDispHeat = .false.   !散逸加熱を考慮するかのスイッチ
     FlagDExnerDtTurb =.true. !圧力方程式に乱流拡散項を考慮するかのスイッチ
  /


== 雲物理パラメタリゼーションの設定

=== FlagCloudMethod = "K1969"

暖かい雨のパラメタリゼーション (Kessler, 1969) を用いる場合の設定. 

  &cloudphys_k1969_nml
    Planet            = ""        ! "Earth" of "Jupiter"
                                  !   Earth:   FactorJ = 1.0 に設定される
                                  !   Jupiter: FactorJ = 3.0 に設定される
    FactorJ           = 1.0d0     ! 雲物理過程のパラメータ
                                  !   木星では 3.0d0
                                  !   地球では 1.0d0 とする
    AutoConvTime      = 1000.0d0  ! 併合成長の時定数 [sec]
    QMixCr            = 1.0d-3    ! 併合成長を生じる臨界混合比 [kg/kg]
    FlagDExnerDtCloud = .true.    !圧力の式に凝結の効果を考慮するか否か
                                  !考慮しない場合は値を .false. にする.
    FlagDExnerDtFall  = .true.    !圧力の式に落下の効果を考慮するか否か
                                  !考慮しない場合は値を .false. にする.
    FactorFallRain    = 1.0d0     !雨の落下の有無
                                  !考慮しない場合は値をゼロにする.
    FactorCloud2Rain  = 1.0d0     !雲から雨への変換の有無
                                  !考慮しない場合は値をゼロにする.
    FactorRain2Gas    = 1.0d0     !雨から蒸気への変換の有無
                                  !考慮しない場合は値をゼロにする.
    FactorCloud2Gas   = 1.0d0     !雲から蒸気への変換の有無 
                                  !考慮しない場合は値をゼロにする.
  /

=== FlagCloudMethod = "MarsCond"

Yamasita et al (投稿準備中) で用いている拡散成長

  &cloudphys_marscond_nml
    DensIce     = 1.565d3  ! 固相の密度 [kg/m^3]
    NumAerosol  = 0.0d0    ! エアロゾルの数密度 [1/kg]
    RadiAerosol = 0.0d0    ! エアロゾルの数密度 [1/kg]
    Kd          = 0.0d0    ! 大気の熱伝導係数 [W/K m]
    SatRatioCr  = 0.0d0    ! 臨界飽和比 
    SatRtWetAdia = 0.0d0   ! 湿潤断熱線の飽和比 
    CO2LatHeat  = 0.0d0    ! 単位質量あたりの凝結熱 [J/kg]
    CDensCr     = 5.0d-5   ! 閾値
  /

== 放射過程

=== FlagRadMethod = "HeatConst"

簡単放射: とある高度を一様加熱・冷却.

  &radiation_simple_nml
    RadHeatRate= 0.0d0,     !一様放射強制の大きさ [K/day]
    HeightUp   = 10.0d3,    !放射強制を与える鉛直領域の上限
    HeightDown = 0.0d3,     !放射強制を与える鉛直領域の下限 
    FlagDExnerDtRad = .true. !圧力方程式で放射による加熱冷却の寄与を考慮するかのスイッチ
  /

=== FlagRadMethod = "HeatVary"

簡単放射: 地表面から HeightDown までは RadHeatRate で冷却. 
HeightUp より上空は加熱率ゼロになるように加熱率を減少させる. 

  &radiation_simple_nml
    RadHeatRate= 0.0d0,      !放射強制の大きさ [K/day]
    HeightUp   = 10.0d3,     !放射強制を与える鉛直領域の上限
    HeightDown = 0.0d3,      !放射強制を与える鉛直領域の下限 
    FlagDExnerDtRad = .true. !圧力方程式で放射による加熱冷却の寄与を考慮するかのスイッチ
  /

=== FlagRadMethod = "HeatBalance"

簡単放射モジュール: とある高度領域を一様冷却・加熱する. 
冷却率は設定ファイルで与える. 加熱率は冷却と釣り合うようにモデル内部で決める. 

  &radiation_heatbalance_nml
    RadCoolRate = 0.0d0      !一様放射加熱率 [K/day]                 
    HeightHeatUp   = 0.0d0   !加熱領域の上端の高度
    HeightHeatDown = 0.0d0   !加熱領域の下端の高度
    HeightCoolUp   = 0.0d0   !冷却領域の上端の高度
    HeightCoolDown = 0.0d0   !冷却領域の下端の高度
    FlagDExnerDtRad = .true. !圧力方程式で放射による加熱冷却の寄与を考慮するかのスイッチ
  /

=== FlagRadMethod = "Baker1998"

Baker et al. (1998) の金星計算での放射過程. 設定する項目は無い. 

=== FlagRadMethod = "Sounding"

  &radiation_sounding_nml
    SoundingFile    = ""       ! サウンディングファイル
    AltCol          = 0        !「高度」の列番号 (サウンディングファイル内)
    SWaveCol        = 0        !「短波放射による加熱率」の列番号 (サウンディングファイル内)   
    LWaveCol        = 0        !「長波放射による加熱率」の列番号 (サウンディングファイル内)   
    FlagDExnerDtRad = .true.   !圧力方程式で放射による加熱冷却の寄与を考慮するかのスイッチ
  /

== 地表面過程
  
=== FlagSurfaceMethod = "Diff"  

下部境界の拡散係数を決めうちする.   

  &surfaceflux_diff_nml
     Kappa = 800.0d0               ! 下部境界での乱流拡散係数
     FlagDExnerDtSurf = .true.     ! Flag for diabatice heating term in pressure equation
  /


===  FlagSurfaceMethod = "Bulk"  

  &surfaceflux_bulk_nml
    FlagConstBulkCoef
                            ! Flag for using constant bulk coefficient
    FlagUseOfBulkCoefInNeutralCond
                            ! Flag for using bulk coefficient in neutral condition
    FlagDExnerDtSurf = .true.  
                            ! Flag for diabatice heating term in pressure equation
    ConstBulkCoef           
                            ! バルク係数一定値. 
                            ! Steady value of bulk coefficient
    VelMinForRi = 1.0d-8    ! リチャード数計算用速度下限値
                            ! Lower limit of velocity for Ri
    SfcRoughLength = 1.0d-2 ! 祖度長さ
                            ! Roughness length
    Vel0 = 0.0d0            ! 下層での水平速度嵩上げ値
                            ! 
    VelBulkCoefMin = 0.0d0  ! $ u $ バルク係数最小値. 
                            ! Minimum value of $ u $ bulk coefficient
    TempBulkCoefMin = 0.0d0 ! $ T $ バルク係数最小値. 
                            ! Minimum value of $ T $ bulk coefficient
    QmixBulkCoefMin = 0.0d0 ! $ q $ バルク係数最小値. 
                            ! Minimum value of $ q $ bulk coefficient
    VelBulkCoefMax = 1.0d2  ! $ u $ バルク係数最大値. 
                            ! Maximum value of $ u $ bulk coefficient
    TempBulkCoefMax = 1.0d2 ! $ T $ バルク係数最大値. 
                            ! Maximum value of $ T $ bulk coefficient
    QmixBulkCoefMax = 1.0d2 ! $ q $ バルク係数最大値. 
                            ! Maximum value of $ q $ bulk coefficient
  /

===  FlagSurfaceMethod = "Const"

上部, 下部境界からの一定の熱・運動量・物質フラックスを与えた場合には
以下の項目を設定する. 

  &surfaceflux_const_nml
     SfcXMomFluxBtm = 0.0d0    ! X 方向の運動量フラックス (下部境界)
     SfcXMomFluxTop = 0.0d0    ! X 方向の運動量フラックス (上部境界)
     SfcYMomFluxBtm = 0.0d0    ! Y 方向の運動量フラックス (下部境界)
     SfcYMomFluxTop = 0.0d0    ! Y 方向の運動量フラックス (上部境界)
     SfcHeatFluxBtm = 0.0d0    ! 熱フラックス (下部境界)
     SfcHeatFluxTop = 0.0d0    ! 熱フラックス (上部境界)
     SfcQmixFluxBtm = 0.0d0    ! 物質フラックス (下部境界)
     SfcQmixFluxTop = 0.0d0    ! 物質フラックス (上部境界)
     FlagDExnerDtSurf = .true. ! Flag for diabatice heating term in pressure equation
  /

===  FlagSurfaceMethod = "Baker1998"

Baker et al. (1998) の金星計算での地表面(?)過程. 設定する項目は無い. 


=end
