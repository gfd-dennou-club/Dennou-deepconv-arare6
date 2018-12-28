#= Makefile for deepconv/arare reference mannual.
#
# Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
# Version::   $Id: code_reference.rd,v 1.2 2011/12/19 08:05:04 odakker Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:../../COPYRIGHT]
#
#
########################################################################
#

=begin JA

= deepconv/arare ソースコードリファレンス

#* 履歴
#  * 2011/12/19(小高正嗣) 更新
#  * 2008/06/19(山下達也) 更新
#  * 2007/10/19(小高正嗣) 更新
#  * 2006/11/10(小高正嗣) 更新
#  * 2006/11/09(小高正嗣) 新規作成

=end JA
=begin EN

= Deepconv/arare source code reference

#* History
#  * 2008/06/19 (Tatsuya Yamashita) Update
#  * 2007/10/19 (Masatsugu Odaka) Update
#  * 2006/11/10 (Masatsugu Odaka) Update
#  * 2006/11/09 (Masatsugu Odaka) Initial release

=end EN


=begin JA

== 実行プログラムの一覧

((<主プログラム(main)|URL:../src/main>)):
* ((<arare.f90|URL:code_reference/files/__/src/main/arare_f90.html>)): 
  主プログラム
* ((<arare_init-data.f90|URL:code_reference/files/__/src/main/arare_init-data_f90.html>)): 
  初期値計算用プログラム



=end JA
=begin EN

== Execute programs

((<Main program (main)|URL:../src/main>)):
* ((<arare.f90|URL:code_reference/files/main/arare_f90.html>)): 
  Main program
* ((<arare_init-data.f90|URL:code_reference/files/__/src/main/arare_init-data_f90.html>)): 
  Main program for initial data calculation.


=end EN


=begin JA
== サブルーチンとモジュールの一覧

((<化学過程(chemdata)|URL:../src/chemdata>))
* ((<ChemData|URL:code_reference/classes/ChemData.html>)):
  化学種データ保管モジュール


((<力学過程(dynamic)|URL:../src/dynamics>))
* ((<DynamicsHEVI|URL:code_reference/classes/DynamicsHEVI.html>)):
  力学過程計算用モジュール


((<初期値設定(env)|URL:../src/env/>))
* ((<initialdata_basic|URL:code_reference/classes/initialdata_basic.html>)):
  基本場設定用モジュール
* ((<initialdata_disturb|URL:code_reference/classes/initialdata_disturb.html>)):
  初期擾乱設定用モジュール

* ((<initialdata_yamasaki1983|URL:code_reference/classes/initialdata_yamasaki1983.html>)):
  Yamasaki (1983) の初期値設定用モジュール
* ((<initialdata_Toon2002|URL:code_reference/classes/initialdata_Toon2002.html>)):
  Toon (2002) の初期値設定用モジュール
* ((<initialdata_takemi2007|URL:code_reference/classes/initialdata_takemi2007.html>)):
  Takemi (2007) の初期値設定用モジュール


((<入出力(io)|URL:../src/io/>))
* ((<ReStartFileIO|URL:code_reference/classes/ReStartFileIO.html>)):
  リスタートファイルへの出力モジュール
* ((<HistoryFileIO|URL:code_reference/classes/HistoryFileIO.html>)):
  ヒストリファイルへの出力モジュール
* ((<Arare4fileio|URL:code_reference/classes/Arare4fileio.html>)):
  deepconv/arare4 で作成したリスタートファイルの入力モジュール


((<物理過程(physics)|URL:../src/physics/>))
* ((<cloudphys_k1969|URL:code_reference/classes/cloudphys_k1969.html>)):
  Kessler (1969) による雲微物理パラメタリゼーション計算用モジュール
* ((<Cloudphys_MarsCond|URL:code_reference/classes/Cloudphys_MarsCond.html>)):
  火星大気雲微物理パラメタリゼーション計算用モジュール
* ((<ECCM|URL:code_reference/classes/ECCM.html>)):
  上昇する断熱気塊の温度減率と, その場合の平衡大気構造を計算するモジュール
* ((<MoistAdjust|URL:code_reference/classes/MoistAdjust.html>)):
  湿潤飽和調節法モジュール
* ((<Radiation_HeatBalance|URL:code_reference/classes/Radiation_HeatBalance.html>)):
  放射を模した熱強制計算用モジュール
* ((<Radiation_Simple|URL:code_reference/classes/Radiation_Simple.html>)):
  放射を模した熱強制計算用モジュール(単純バージョン)
* ((<Surfaceflux_bulk|URL:code_reference/classes/Surfaceflux_bulk.html>)):
  バルク法による地表フラックス計算用モジュール
* ((<Surfaceflux_diff|URL:code_reference/classes/Surfaceflux_diff.html>)):
  拡散による地表フラックス計算用モジュール
* ((<Turbulence_kw1978|URL:code_reference/classes/Turbulence_kw1978.html>)):
  Klemp & Wilhelmson (1978) の乱流パラメタリゼーション計算用モジュール


((<初期設定(setup)|URL:../src/setup/>))
* ((<argset|URL:code_reference/classes/argset.html>)):
  コマンドライン引数解釈用モジュール
* ((<axesset|URL:code_reference/classes/axesset.html>)):
  3 次元等間隔交互格子格子点設定モジュール
* ((<basicset|URL:code_reference/classes/basicset.html>)):
  基本場設定モジュール
* ((<ChemCalc|URL:code_reference/classes/ChemCalc.html>)):
  化学関連量計算モジュール
* ((<clockset|URL:code_reference/classes/clockset.html>)):
  計算時間情報取得用モジュール
* ((<composition|URL:code_reference/classes/composition.html>)):
  化学定数設定モジュール
* ((<constats|URL:code_reference/classes/constants.html>)):
  定数設定用モジュール  
* ((<constats0|URL:code_reference/classes/constants0.html>)):
  物理・数学定数設定モジュール
* ((<dataset|URL:code_reference/classes/dataset.html>)):
  物理・化学的パラメータ設定モジュール 
* ((<fileset|URL:code_reference/classes/fileset.html>)):
  入出力ファイル名設定モジュール
* ((<gridset|URL:code_reference/classes/gridset.html>)):
  格子点配列サイズ設定モジュール
* ((<mpi_wrapper|URL:code_reference/classes/mpi_wrapper.html>)):
  MPI ラッパーモジュール
* ((<namelist_util|URL:code_reference/classes/namelist_util.html>)):
  NAMELIST ファイル入力に関するモジュール
* ((<timset|URL:code_reference/classes/timeset.html>)):
  時間積分用パラメータ設定モジュール


((<下請けモジュール(util)|URL:../src/util/>))
* ((<CFLCheck|URL:code_reference/classes/CFLCheck.html>)):
  CFL 条件確認モジュール
* ((<Damping|URL:code_reference/classes/Damping.html>)):
  音波減衰項とスポンジ層での摩擦項の計算モジュール
* ((<FillNegative|URL:code_reference/classes/FillNegative.html>)):
  雲水量などの正定値量の穴埋め計算モジュール
* ((<setmargin|URL:code_reference/classes/setmargin.html>)):
  糊代領域配列サイズ設定モジュール
* ((<TimeFilter|URL:code_reference/classes/TimeFilter.html>)):
  時間フィルター計算モジュール
* ((<xyz_bc_module|URL:code_reference/classes/xyz_bc_module.html>)):
  境界条件設定モジュール
* ((<xyz_deriv_c4_module|URL:code_reference/classes/xyz_deriv_c4_module.html>))
  4 次精度中心差分計算モジュール
* ((<xys_deriv_module|URL:code_reference/classes/xyz_deriv_module.html>)):
  2 次精度中心差分計算モジュール


=end JA
=begin EN

== Subroutines and modules

((<Chemical process (chemdata)|URL:../src/chemdata>))
* ((<ChemData|URL:code_reference/classes/ChemData.html>)):
  Chemical data module


((<Dynamics (dynamic)|URL:../src/dynamics>))
* ((<DynamicsHEVI|URL:code_reference/classes/DynamicsHEVI.html>)):
  Module for dynamical process

((<Initial environment setup (env)|URL:../src/env/>))
* ((<initialdata_basic|URL:code_reference/classes/initialdata_basic.html>)):
  Basic state set up module
* ((<initialdata_disturb|URL:code_reference/classes/initialdata_disturb.html>)):
  Initial disturbance set up module

* ((<initialdata_yamasaki1983|URL:code_reference/classes/initialdata_yamasaki1983.html>)):
  Initial value used by Yamasaki (1983) setup module
* ((<initialdata_Toon2002|URL:code_reference/classes/initialdata_Toon2002.html>)):
  Initial value used by Toon et al. (2002) setup module
* ((<initialdata_takemi2007|URL:code_reference/classes/initialdata_takemi2007.html>)):
  Initial value used by Takemi (2007) setup module


((<"Input/Output (io)"|URL:../src/io/>))
* ((<HistoryFileIO|URL:code_reference/classes/HistoryFileIO.html>)):
  I/O module of history files
* ((<ReStartFileIO|URL:code_reference/classes/ReStartFileIO.html>)):
  I/O module of restart file
* ((<Arare4fileio|URL:code_reference/classes/Arare4fileio.html>)):
  I/O module of restart file generated by deepconv/arare4


((<Physics (physics)|URL:../src/physics/>))
* ((<cloudphys_k1969|URL:code_reference/classes/cloudphys_k1969.html>)):
  Kessler (1969) cloud parameterization module
* ((<Cloudphys_MarsCond|URL:code_reference/classes/Cloudphys_MarsCond.html>)):
  Cloud parameterization module for the Martian atmosphere
* ((<ECCM|URL:code_reference/classes/ECCM.html>)):
  ECCM (Ensemble Cloud Condensation Model) module
* ((<MoistAdjust|URL:code_reference/classes/MoistAdjust.html>)):
  Moist adjustment module
* ((<Radiation_HeatBalance|URL:code_reference/classes/Radiation_HeatBalance.html>)):
  Thermal forcing module associated with atmospheric radiation
* ((<Radiation_Simple|URL:code_reference/classes/Radiation_Simple.html>)):
  Thermal forcing module associated with atmospheric radiation (simple verion)
* ((<Surfaceflux_bulk|URL:code_reference/classes/Surfaceflux_bulk.html>)):
  Surface flux calculation module by using bulk method
* ((<Surfaceflux_diff|URL:code_reference/classes/Surfaceflux_diff.html>)):
  Surface flux calculation module by using diffusion
* ((<Turbulence_kw1978|URL:code_reference/classes/Turbulence_kw1978.html>)):
  Klemp & Wilhelmson (1978) turbulent parameterization module

((<Set up (setup)|URL:../src/setup/>))
* ((<argset|URL:code_reference/classes/argset.html>)):
  Command line argument module
* ((<axesset|URL:code_reference/classes/axesset.html>)):
  Grid arrangement set up module
* ((<basicset|URL:code_reference/classes/basicset.html>)):
  Basic state set up module
* ((<ChemCalc|URL:code_reference/classes/ChemCalc.html>)):
  Chemical process module
* ((<clockset|URL:code_reference/classes/clockset.html>)):
  Clock information set up module
* ((<composition|URL:code_reference/classes/composition.html>)):
  Chemical constants set up module
* ((<constats|URL:code_reference/classes/constants.html>)):
  Constant values set up module
* ((<constats0|URL:code_reference/classes/constants0.html>)):
  Physical and mathmatical constants set up module
* ((<dataset|URL:code_reference/classes/dataset.html>)):
  Physical and chemical parameter set up module
* ((<fileset|URL:code_reference/classes/fileset.html>)):
  I/O file names set up module
* ((<gridset|URL:code_reference/classes/gridset.html>)):
  Grid array size set up module
* ((<mpi_wrapper|URL:code_reference/classes/mpi_wrapper.html>)):
  MPI wrappaer module
* ((<namelist_util|URL:code_reference/classes/namelist_util.html>)):
  NAMELIST file name parameter set up module
* ((<timset|URL:code_reference/classes/timeset.html>)):
  Time integration parameters set up module


((<Utility (util)|URL:../src/util/>))
* ((<CFLCheck|URL:code_reference/classes/CFLCheck.html>)):
  CFL condition check 
* ((<Damping|URL:code_reference/classes/Damping.html>)):
  Sound wave damping term and Rayleigh damping near the upper boundary
* ((<FillNegative|URL:code_reference/classes/FillNegative.html>)):
  Fulfill negative value of positive definite variables
* ((<setmargin|URL:code_reference/classes/setmargin.html>)):
  Array size of margine area set up module
* ((<TimeFilter|URL:code_reference/classes/TimeFilter.html>)):
  Time filter for time integration
* ((<xyz_bc_module|URL:code_reference/classes/xyz_bc_module.html>)):
  Adapting boundary condition 
* ((<xyz_deriv_c4_module|URL:code_reference/classes/xyz_deriv_c4_module.html>))
  4th order centered differentiate scheme
* ((<xys_deriv_module|URL:code_reference/classes/xyz_deriv_module.html>)):
  2'nd order centered differentiate scheme 

=end EN
