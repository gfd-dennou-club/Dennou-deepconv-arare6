!= Program ArareInitData
!
! Authors::   SUGIYAMA Ko-ichiro (杉山耕一朗), ODAKA Masatsugu (小高正嗣)
! Version::   $Id: arare_init-data.f90,v 1.30 2014/07/08 01:01:45 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

program ArareInitData
  !
  ! 非静力学モデル deepconv/arare
  ! 初期値作成用プログラム
  !

  !----- モジュール読み込み ------

  !-----   型宣言, 文字列処理   ----
  use dc_types,       only : STRING, DP
  use dc_string,      only : StoA
  use gtool_history,  only : HistoryPut

  !-----   メッセージ出力   -----
  use dc_message,     only: MessageNotify, MessageSuppressMPI

  !  コマンドライン引数解釈
  use argset,        only : argset_init  

  !  MPI の初期化
  !  Initialize MPI wrapper
  !
  use mpi_wrapper,   only : MPIWrapperInit, MPIWrapperFinalize

  !-----    管理モジュール   -----
  !  入出力ファイル名管理モジュール
  use fileset,    only : fileset_init

  !  格子点管理モジュール 
  use gridset,    only : gridset_init, &
    &                    imin, imax, jmin, jmax, kmin, kmax, &
    &                    nx, ny, nz,  ncmax

  !  基本場設定モジュール
  use basicset,      only : basicset_init

  use axesset,       only : axesset_init

  use constants,     only : constants_init, molwtdry, &
    &                       Grav, TempSfc, PressSfc, &
    &                       pressbasis, gasrdry, cpdry, cvdry

  !  湿潤ルーチン設定モジュール
  use composition,      only: composition_init, molwtwet, SpcWetMolFr, &
    &                         CondNum, IdxCC, IdxCG, SpcWetID

  ! 化学計算
  use ChemCalc, only: ChemCalc_init, xyz_SvapPress

  !  境界条件適用モジュール
  use setmargin

  use axesset, only:  z_Z,          &!スカラー格子点での高度
    &                 dz             !Z 方向の格子点間隔

  use namelist_util, only: NmlutilInit, namelist_filename

  use eccm,     only: ECCM_Dry,      &!
    &                 ECCM_Wet              

  use initialdata_disturb, only:                                  &
    &  initialdata_disturb_random,                                &
    &  initialdata_disturb_circleXZ,                              &
    &  initialdata_disturb_gaussXZ, initialdata_disturb_gaussXY,  &
    &  initialdata_disturb_gaussYZ, initialdata_disturb_gaussXYZ, &
    &  initialdata_disturb_cosXZ,   initialdata_disturb_cosXY,    &
    &  initialdata_disturb_cosYZ,   initialdata_disturb_cosXYZ,   &
    &  initialdata_disturb_coneXZ,  initialdata_disturb_coneXY,   &
    &  initialdata_disturb_coneYZ,                                &
    &  initialdata_disturb_tanh,    initialdata_disturb_tanh_sin, &
    &  initialdata_disturb_dryreg,  initialdata_disturb_moist,    &
    &  initialdata_disturb_square
  use initialdata_baker1998, only: initialdata_baker1998_basic
  use initialdata_Skamarock1994, only: initialdata_SK1994_basic, &
    &                                  initialdata_SK1994_disturbXZ, &
    &                                  initialdata_SK1994_disturbYZ
  use initialdata_yamasaki1983, only: initialdata_yamasaki1983_basic
  use initialdata_takemi2007, only: initialdata_takemi2007_init, &
    &  initialdata_takemi2007_basic, initialdata_takemi2007_wind
  use initialdata_toon2002, only: initialdata_toon2002_basic
  use initialdata_sounding, only: initialdata_sounding_init, &
    & initialdata_sounding_basic, initialdata_sounding_wind

  !-----    入出力モジュール   -----
  !  リスタートファイル入出力モジュール
  use RestartFileIO,    only : ReStartFileIO_Init, ReStartFileIO_Finalize, rstat
  use Arare4InitFileIO, only : Arare4InitFileIO_Init, Arare4InitFileIO_Var_Get, Arare4InitFileIO_BZ_Get, &
    &                          Arare4InitFileIO_MMC_Var_Get, Arare4InitFileIO_MMC_BZ_Get

  !暗黙の型宣言禁止
  implicit none

  !内部変数
  real(DP), allocatable :: pyz_VelX(:,:,:)
  real(DP), allocatable :: xqz_VelY(:,:,:)
  real(DP), allocatable :: xyr_VelZ(:,:,:)
  real(DP), allocatable :: xyz_Exner(:,:,:)
  real(DP), allocatable :: xyz_PTemp(:,:,:)
  real(DP), allocatable :: xyz_Km(:,:,:)
  real(DP), allocatable :: xyz_Kh(:,:,:)
  real(DP), allocatable :: xyz_CDens(:,:,:)
  real(DP), allocatable :: xyzf_QMix(:,:,:,:)
  real(DP), allocatable :: xyz_DensBZ(:,:,:)
  real(DP), allocatable :: xyz_PressBZ(:,:,:)
  real(DP), allocatable :: xyz_ExnerBZ(:,:,:)
  real(DP), allocatable :: xyz_TempBZ(:,:,:)
  real(DP), allocatable :: xyz_PTempBZ(:,:,:)
  real(DP), allocatable :: xyz_VPTempBZ(:,:,:)
  real(DP), allocatable :: xyz_VelSoundBZ(:,:,:)
  real(DP), allocatable :: xyzf_QMixBZ(:,:,:,:)
  real(DP), allocatable :: xyz_EffMolWtBZ(:,:,:)
  real(DP), allocatable :: z_TempBZ(:)
  real(DP), allocatable :: z_PressBZ(:)
  real(DP), allocatable :: xyzf_MolFr(:,:,:,:)
  real(DP), allocatable :: zf_MolFr(:,:)
  real(DP), allocatable :: xyzf_QMixDivMolWt(:,:,:,:)
  real(DP), allocatable :: xyzf_HumBZ(:,:,:,:)
  real(DP)              :: Time
  integer               :: i, j, k, s
  logical, parameter    :: FlagInitData = .true.
  
  !変数定義
  ! Moist 用
  real(DP), save :: Humidity = 0.0d0 !相対湿度
  ! Gauss 用, cos 用, cone 用
  real(DP), save :: PTempMax = 0.0d0 !最大値
  real(DP), save :: ExnerMax = 0.0d0 !最大値
  real(DP), save :: QMixMax = 0.0d0  !最大値
  real(DP), save :: Xc = 0.0d0       !擾乱の中心位置(X方向)
  real(DP), save :: Yc = 0.0d0       !擾乱の中心位置(Y方向)
  real(DP), save :: Zc = 0.0d0       !擾乱の中心位置(鉛直方向)
  real(DP), save :: Xr = 0.0d0       !擾乱の半径(X方向)
  real(DP), save :: Yr = 0.0d0       !擾乱の半径(Y方向)
  real(DP), save :: Zr = 0.0d0       !擾乱の半径(鉛直方向)
  ! tanh 用
  real(DP), save :: PTempMean= 0.0d0 !シアーの中心の値
  real(DP), save :: VelMean  = 0.0d0 !シアーの中心の値
  real(DP), save :: PTempDel = 0.0d0 !シアーでの変化量
  ! Therma-Random 用
  real(DP), save :: Zpos = 0.0d0     !擾乱の Z 座標 [m] 
  ! Square 用
  real(DP), save :: XposMin = 0.0d0    !乾燥域の X 座標 [m] 
  real(DP), save :: YposMin = 0.0d0    !乾燥域の Y 座標 [m] 
  real(DP), save :: ZposMin = 0.0d0    !乾燥域の Z 座標 [m] 
  real(DP), save :: XposMax = 0.0d0    !乾燥域の X 座標 [m] 
  real(DP), save :: YposMax = 0.0d0    !乾燥域の Y 座標 [m] 
  real(DP), save :: ZposMax = 0.0d0    !乾燥域の Z 座標 [m] 
  ! WindConst 用
  real(DP), save :: VelX0 = 0.0d0   !初期に一定の速度を与える場合
  real(DP), save :: VelY0 = 0.0d0   !初期に一定の速度を与える場合
  real(DP), save :: VelZ0 = 0.0d0   !初期に一定の速度を与える場合

  real(DP), save :: TempTr  = 10.0d0    !対流圏界面の温度 [k]
  real(DP), save :: Dheight = 10.0d3    !重み関数のパラメータ [m]
!  real(DP), save :: HeightTr= 10000.0d3 !対流圏界面の高度 [m]

  integer               :: IDBasic                = 0
  integer, parameter    :: IDBasicArare4          = 1
  integer, parameter    :: IDBasicArare4mmc       = 2
  integer, parameter    :: IDBasicYamasaki1983    = 3
  integer, parameter    :: IDBasicTakemi2007      = 4
  integer, parameter    :: IDBasicIsoThermal      = 5
  integer, parameter    :: IDBasicDry             = 6
  integer, parameter    :: IDBasicMoist           = 7
  integer, parameter    :: IDBasicToon2002        = 8
  integer, parameter    :: IDBasicSounding        = 9
  integer, parameter    :: IDBasicBaker1998       = 10
  integer, parameter    :: IDBasicSK1994   = 11
  integer               :: IDDisturbPTemp         = 0
  integer, parameter    :: IDDisturbPTempGaussXY  = 1
  integer, parameter    :: IDDisturbPTempGaussXZ  = 2
  integer, parameter    :: IDDisturbPTempGaussYZ  = 3
  integer, parameter    :: IDDisturbPTempGaussXYZ = 4
  integer, parameter    :: IDDisturbPTempRandom   = 5
  integer, parameter    :: IDDisturbPTempSquare= 6
  integer, parameter    :: IDDisturbPTempCosXY    = 7
  integer, parameter    :: IDDisturbPTempCosXZ    = 8
  integer, parameter    :: IDDisturbPTempCosYZ    = 9
  integer, parameter    :: IDDisturbPTempCosXYZ   = 10
  integer, parameter    :: IDDisturbPTempConeXY   = 11
  integer, parameter    :: IDDisturbPTempConeXZ   = 12
  integer, parameter    :: IDDisturbPTempConeYZ   = 13
  integer, parameter    :: IDDisturbPTempTanh     = 14
  integer, parameter    :: IDDisturbPTempSK1994XZ = 15
  integer, parameter    :: IDDisturbPTempSK1994YZ = 16
  integer, parameter    :: IDDisturbPTempCircleXZ = 17
  integer               :: IDDisturbExner         = 0
  integer, parameter    :: IDDisturbExnerGaussXY  = 1
  integer, parameter    :: IDDisturbExnerGaussXZ  = 2
  integer, parameter    :: IDDisturbExnerGaussYZ  = 3
  integer, parameter    :: IDDisturbExnerGaussXYZ = 4
  integer, parameter    :: IDDisturbExnerCosXY    = 5
  integer, parameter    :: IDDisturbExnerCosXZ    = 6
  integer, parameter    :: IDDisturbExnerCosYZ    = 7
  integer, parameter    :: IDDisturbExnerCosXYZ   = 8
  integer, parameter    :: IDDisturbExnerConeXY   = 9
  integer, parameter    :: IDDisturbExnerConeXZ   = 10
  integer, parameter    :: IDDisturbExnerConeYZ   = 11
  integer               :: IDDisturbQMix          = 0
  integer, parameter    :: IDDisturbQMixGaussXY   = 1
  integer, parameter    :: IDDisturbQMixGaussXZ   = 2
  integer, parameter    :: IDDisturbQMixGaussYZ   = 3
  integer, parameter    :: IDDisturbQMixGaussXYZ  = 4
  integer, parameter    :: IDDisturbQMixDryreg    = 5
  integer, parameter    :: IDDisturbQMixMoist     = 6
  integer, parameter    :: IDDisturbQMixCosXY     = 7
  integer, parameter    :: IDDisturbQMixCosXZ     = 8
  integer, parameter    :: IDDisturbQMixCosYZ     = 9
  integer, parameter    :: IDDisturbQMixCosXYZ    = 10
  integer, parameter    :: IDDisturbQMixConeXY    = 11
  integer, parameter    :: IDDisturbQMixConeXZ    = 12
  integer, parameter    :: IDDisturbQMixConeYZ    = 13
  integer, parameter    :: IDDisturbQMixSquare = 14
  integer               :: IDDisturbWind          = 0
  integer, parameter    :: IDDisturbWindTakemi2007= 1
  integer, parameter    :: IDDisturbWindSounding  = 2
  integer, parameter    :: IDDisturbWindZonal     = 3
  integer, parameter    :: IDDisturbWindTanh      = 4
  integer               :: IDDisturb              = 0
  integer, parameter    :: IDDisturbArare4        = 1
  integer, parameter    :: IDDisturbArare4mmc     = 2
  

  !------------------------------------------
  ! 初期化手続き ; Initialize procedure 
  !
  call MainInit

  Time = 0.0d0
  call MessageNotify( "M", "main", "Making Initial data...." )

  !---------------------------------------------------------------  
  ! 基本場を作成する. 
  !
  select case(IDBasic)
    
  case (IDBasicArare4)
    ! deepconv/arare4 のヒストリーファイルから初期値を読み込む
    !
    call MessageNotify( "M", "main", "Making Initial data using arare4 files ...." )
    call Arare4Initfileio_BZ_Get( &
      & xyz_PressBZ,   & ! (out)
      & xyz_ExnerBZ,   & ! (out)
      & xyz_TempBZ,    & ! (out)
      & xyz_PTempBZ,   & ! (out)
      & xyz_DensBZ,    & ! (out)
      & xyz_VelSoundBZ,& ! (out)
      & xyzf_QMixBZ,   & ! (out)
      & xyz_EffMolWtBZ & ! (out)
      & )

  case (IDBasicArare4mmc)
    ! deepconv/arare4 火星版のヒストリーファイルから初期値を読み込む
    !    
    call MessageNotify( "M", "main", "Making Initial data using arare4 files ...." )

    call Arare4Initfileio_MMC_BZ_Get( &
      & xyz_PressBZ,   & ! (out)
      & xyz_ExnerBZ,   & ! (out)
      & xyz_TempBZ,    & ! (out)
      & xyz_PTempBZ,   & ! (out)
      & xyz_DensBZ,    & ! (out)
      & xyz_VelSoundBZ,& ! (out)
      & xyzf_QMixBZ,   & ! (out)
      & xyz_EffMolWtBZ & ! (out)
      & )
    
  case(IDBasicIsoThermal) 
    ! 等温/重力考慮しない場合
    !  
    call MessageNotify( "M", "main", "Making Initial data (basic) named IsoThermal..." )
    
    z_TempBZ = TempSfc
    z_PressBZ = PressSfc
    zf_MolFr = 0.0d0
    
    ! 残りの基本場の値を決める
    call DetermineBZ            

  case(IDBasicDry) 
    ! 乾燥断熱的な初期場
    !
    call MessageNotify( "M", "main", "Making Initial data (basic) named DRY..." )
    call eccm_dry( SpcWetMolFr(1:ncmax), Humidity, z_TempBZ, z_PressBZ, zf_MolFr )
    if (minval(z_TempBZ(1:nz)) < TempTr) then 
!      call initialdata_basic_strat( z_TempBZ, z_PressBZ ) !(inout)
      call initialdata_basic_strat_v2( z_TempBZ, z_PressBZ ) !(inout)
    end if
    
    ! 残りの基本場の値を決める
    call DetermineBZ            
    
! デバッグ済んでいない
!   
!  case(IDBasicMoist) 
!    ! 湿潤断熱的な初期場
!    !
!    call MessageNotify( "M", "main", "Making Initial data (basic) named MOIST..." )
!    call eccm_wet( SpcWetMolFr(1:ncmax), Humidity, z_TempBZ, z_PressBZ, zf_MolFr )
!    call initialdata_basic_strat_v2( z_TempBZ, z_PressBZ ) !(inout)
!
!    ! 残りの基本場の値を決める
!    call DetermineBZ            

  case (IDBasicSounding)
    ! サウンディングファイルから値を読み込む
    !
    call MessageNotify( "M", "main", "Making Initial data using sounding files ...." )
    call initialdata_sounding_basic( z_TempBZ, z_PressBZ )  ! (out)
    zf_MolFr = 0.0d0
    
    ! 残りの基本場の値を決める
    call DetermineBZ            

  case(IDBasicToon2002) 
    ! 火星用
    !
    call MessageNotify( "M", "main", "Making Initial data (basic) named Toon et al. 2002..." )
    call initialdata_toon2002_basic( z_TempBZ, z_PressBZ ) !(out)
!    call initialdata_basic_strat_v2( z_TempBZ, z_PressBZ ) !(inout)

    ! 残りの基本場の値を決める
    call DetermineBZ            
    
  case(IDBasicTakemi2007)
    ! Takemi(2007)の基本場を使用する場合
    !
    call MessageNotify( "M", "main", "Making Initial data named Takemi2007..." )
    call initialdata_takemi2007_basic( z_TempBZ, z_PressBZ, zf_MolFr )

    ! 残りの基本場の値を決める
    call DetermineBZ            

  case(IDBasicYamasaki1983) 
    ! Yamasaki(1983)の温度と相対湿度の観測値を使用する場合 
    !
    call MessageNotify( "M", "main", "Making Initial data named Yamasaki1983..." )
    call initialdata_yamasaki1983_basic( z_TempBZ, z_PressBZ, zf_MolFr )
    
    ! 残りの基本場の値を決める
    call DetermineBZ            

  case(IDBasicBaker1998) 
    ! Baker and Schubert (1998) の温度と相対湿度の観測値を使用する場合 
    !
    call MessageNotify( "M", "main", "Making Initial data named Baker and Shubert 1998..." )
    call initialdata_baker1998_basic( z_TempBZ, z_PressBZ ) !(out)
    zf_MolFr = 0.0d0
    
    ! 残りの基本場の値を決める
    call DetermineBZ            

  case(IDBasicSK1994) 
    ! 温位:  Skamarock & Klemp (1994) 
    !
    call MessageNotify( "M", "main", "Making Initial data named Skamarock & Klemp (1994) ..." )
    call initialdata_SK1994_basic( z_TempBZ, z_PressBZ ) !(out)
    zf_MolFr = 0.0d0
    
    ! 残りの基本場の値を決める
    call DetermineBZ            
    
  end select

  !---------------------------------------------------------
  ! 擾乱場の値を決める
  ! 

  ! 温位: 擾乱場の値を決める
  !    
  select case(IDDisturbPTemp)

  case(IDDisturbPTempGaussXY)
    ! 温位: ガウシアンな分布を与える (XY)
    !          
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXY..." )
    call initialdata_disturb_gaussXY(PTempMax, Xc, Xr, Yc, Yr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempGaussXZ)
    ! 温位: ガウシアンな分布を与える (XZ)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXZ..." )
    call initialdata_disturb_gaussXZ(PTempMax, Xc, Xr, Zc, Zr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempGaussYZ)
    ! 温位: ガウシアンな分布を与える (YZ)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussYZ..." )
    call initialdata_disturb_gaussYZ(PTempMax, Yc, Yr, Zc, Zr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempGaussXYZ)
    ! 温位: ガウシアンな分布を与える (XYZ)
    !                   
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXYZ..." )
    call initialdata_disturb_gaussXYZ(PTempMax, Xc, Xr, Yc, Yr, Zc, Zr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempRandom)
    ! 温位: 指定された高度にランダムな分布を与える. 
    !                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named random..." )
    call initialdata_disturb_random(PTempMax, Zpos, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempSquare)
    ! 温位: 矩形な擾乱を与える
    !                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named Square..." )
    call initialdata_disturb_square(PTempMax, XposMin, XposMax, YposMin, YposMax, ZposMin, ZposMax, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempCosXY)
    ! 温位: A [cos(πL) + 1]*0.5 ( L < 1.0 ) or 0.0な分布を与える (XY)
    !          
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXY..." )
    call initialdata_disturb_cosXY(PTempMax, Xc, Xr, Yc, Yr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempCosXZ)
    ! 温位: A [cos(πL) + 1]*0.5 ( L < 1.0 ) or 0.0な分布を与える (XZ)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXZ..." )
    call initialdata_disturb_cosXZ(PTempMax, Xc, Xr, Zc, Zr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempCosYZ)
    ! 温位: A [cos(πL) + 1]*0.5 ( L < 1.0 ) or 0.0な分布を与える (YZ)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosYZ..." )
    call initialdata_disturb_cosYZ(PTempMax, Yc, Yr, Zc, Zr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempCosXYZ)
    ! 温位: A [cos(πL) + 1]*0.5 ( L < 1.0 ) or 0.0な分布を与える (XYZ)
    !                   
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXYZ..." )
    call initialdata_disturb_cosXYZ(PTempMax, Xc, Xr, Yc, Yr, Zc, Zr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )
    
  case(IDDisturbPTempConeXY)
    ! 温位: 円錐な分布を与える (XY)
    !          
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeXY..." )
    call initialdata_disturb_coneXY(PTempMax, Xc, Xr, Yc, Yr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempConeXZ)
    ! 温位: 円錐な分布を与える (XZ)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeXZ..." )
    call initialdata_disturb_coneXZ(PTempMax, Xc, Xr, Zc, Zr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempConeYZ)
    ! 温位: 円錐な分布を与える (YZ)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeYZ..." )
    call initialdata_disturb_coneYZ(PTempMax, Yc, Yr, Zc, Zr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempTanh)
    ! 温位: tanh な分布を与える (YZ)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named tanh..." )
    call initialdata_disturb_tanh_sin(PTempMean, PTempDel, Zc, Zr, xyz_PTemp, xyz_PTempBZ)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempSK1994XZ)
    ! 温位:  Skamarock & Klemp (1994)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named SK1994..." )
    call initialdata_SK1994_disturbXZ(PTempMax, Xc, Xr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempSK1994YZ)
    ! 温位:  Skamarock & Klemp (1994)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named SK1994..." )
    call initialdata_SK1994_disturbYZ(PTempMax, Yc, Yr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempCircleXZ)
    ! 温位: 円領域を同じ温位擾乱 (XZ)
    !          
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CircleXZ..." )
    call initialdata_disturb_circleXZ(PTempMax, Xc, Xr, Zc, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  end select

  
  ! エクスナー関数: 擾乱場の値を決める    
  !                
  select case(IDDisturbExner)

  case(IDDisturbExnerGaussXY)
    ! エクスナー関数: ガウシアンな分布を与える (XY)
    !                         
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXY..." )
    call initialdata_disturb_gaussXY(ExnerMax, Xc, Xr, Yc, Yr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerGaussXZ)
    ! エクスナー関数: ガウシアンな分布を与える (XZ)
    !                            
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXZ..." )
    call initialdata_disturb_gaussXZ(ExnerMax, Xc, Xr, Zc, Zr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerGaussYZ)
    ! エクスナー関数: ガウシアンな分布を与える (YZ)
    !                            
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussYZ..." )
    call initialdata_disturb_gaussYZ(ExnerMax, Yc, Yr, Zc, Zr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerGaussXYZ)      
    ! エクスナー関数: ガウシアンな分布を与える (XYZ)
    !                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXYZ..." )
    call initialdata_disturb_gaussXYZ(ExnerMax, Xc, Xr, Yc, Yr, Zc, Zr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerCosXY)
    ! エクスナー関数: A [cos(πL) + 1]*0.5 ( L < 1.0 ) or 0.0な分布を与える (XY)
    !                         
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXY..." )
    call initialdata_disturb_cosXY(ExnerMax, Xc, Xr, Yc, Yr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerCosXZ)
    ! エクスナー関数: A [cos(πL) + 1]*0.5 ( L < 1.0 ) or 0.0な分布を与える (XZ)
    !                            
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXZ..." )
    call initialdata_disturb_cosXZ(ExnerMax, Xc, Xr, Zc, Zr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerCosYZ)
    ! エクスナー関数: A [cos(πL) + 1]*0.5 ( L < 1.0 ) or 0.0な分布を与える (YZ)
    !                            
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosYZ..." )
    call initialdata_disturb_cosYZ(ExnerMax, Yc, Yr, Zc, Zr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerCosXYZ)      
    ! エクスナー関数: A [cos(πL) + 1]*0.5 ( L < 1.0 ) or 0.0な分布を与える (XYZ)
    !                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXYZ..." )
    call initialdata_disturb_cosXYZ(ExnerMax, Xc, Xr, Yc, Yr, Zc, Zr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerConeXY)
    ! エクスナー関数: 円錐な分布を与える (XY)
    !                         
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeXY..." )
    call initialdata_disturb_coneXY(ExnerMax, Xc, Xr, Yc, Yr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerConeXZ)
    ! エクスナー関数: 円錐な分布を与える (XZ)
    !                            
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeXZ..." )
    call initialdata_disturb_coneXZ(ExnerMax, Xc, Xr, Zc, Zr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerConeYZ)
    ! エクスナー関数: 円錐な分布を与える (YZ)
    !                            
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeYZ..." )
    call initialdata_disturb_coneYZ(ExnerMax, Yc, Yr, Zc, Zr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  end select


  ! 混合比: 擾乱場の値を決める
  !                    
  select case(IDDisturbQMix)

  case(IDDisturbQMixGaussXY)
    ! 混合比: ガウシアンな分布を与える (XY)
    !                               
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXY..." )
    call initialdata_disturb_gaussXY(QMixMax, Xc, Xr, Yc, Yr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )
    
  case(IDDisturbQMixGaussXZ)
    ! 混合比: ガウシアンな分布を与える (XZ)
    !                                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXZ..." )
    call initialdata_disturb_gaussXZ(QMixMax, Xc, Xr, Zc, Zr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixGaussYZ)
    ! 混合比: ガウシアンな分布を与える (YZ)
    !                                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussYZ..." )
    call initialdata_disturb_gaussYZ(QMixMax, Yc, Yr, Zc, Zr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixGaussXYZ)
    ! 混合比: ガウシアンな分布を与える (XYZ)
    !                                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXYZ..." )
    call initialdata_disturb_gaussXYZ(QMixMax, Xc, Xr, Yc, Yr, Zc, Zr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixDryreg)
    ! 混合比: 乾燥領域を作る
    !                            
    call MessageNotify( "M", "main", "Making Initial data (disturb) named Dryreg..." )
    call initialdata_disturb_Dryreg(XposMin, XposMax, YposMin, YposMax, ZposMin, ZposMax, xyzf_QMix)
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixMoist)
    ! 混合比: 湿潤断熱的な分布を与える 
    !                                  
    call MessageNotify( "M", "main", "Making Initial data (disturb) named MOIST..." )
    call initialdata_disturb_moist(Humidity, xyzf_QMix)
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixCosXY)
    ! 混合比: A [cos(πL) + 1]*0.5 ( L < 1.0 ) or 0.0な分布を与える (XY)
    !                               
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXY..." )
    call initialdata_disturb_cosXY(QMixMax, Xc, Xr, Yc, Yr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )
    
  case(IDDisturbQMixCosXZ)
    ! 混合比: A [cos(πL) + 1]*0.5 ( L < 1.0 ) or 0.0な分布を与える (XZ)
    !                                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXZ..." )
    call initialdata_disturb_cosXZ(QMixMax, Xc, Xr, Zc, Zr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixCosYZ)
    ! 混合比: A [cos(πL) + 1]*0.5 ( L < 1.0 ) or 0.0な分布を与える (YZ)
    !                                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosYZ..." )
    call initialdata_disturb_cosYZ(QMixMax, Yc, Yr, Zc, Zr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixCosXYZ)
    ! 混合比: A [cos(πL) + 1]*0.5 ( L < 1.0 ) or 0.0な分布を与える (XYZ)
    !                                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXYZ..." )
    call initialdata_disturb_cosXYZ(QMixMax, Xc, Xr, Yc, Yr, Zc, Zr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )
    
  case(IDDisturbQMixConeXY)
    ! 混合比: 円錐な分布を与える (XY)
    !                               
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeXY..." )
    call initialdata_disturb_coneXY(QMixMax, Xc, Xr, Yc, Yr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )
    
  case(IDDisturbQMixConeXZ)
    ! 混合比: 円錐な分布を与える (XZ)
    !                                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeXZ..." )
    call initialdata_disturb_coneXZ(QMixMax, Xc, Xr, Zc, Zr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixConeYZ)
    ! 混合比: 円錐な分布を与える (YZ)
    !                                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeYZ..." )
    call initialdata_disturb_coneYZ(QMixMax, Yc, Yr, Zc, Zr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixSquare)
    ! 混合比: 矩形な擾乱を与える
    !                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named Square..." )
    call initialdata_disturb_square(QMixMax, XposMin, XposMax, YposMin, YposMax, ZposMin, ZposMax, xyzf_QMix)
    call SetMargin_xyzf( xyzf_QMix )

  end select

  ! 速度場の値を決める
  !                    
  select case(IDDisturbWind)  
  case (IDDisturbWindTakemi2007)  
    ! Takemi (2007) の風速分布
    !
    call MessageNotify( "M", "main", "Making Initial wind data (disturb) named takemi2007..." )
    call initialdata_takemi2007_wind( pyz_VelX ) 
    call SetMargin_pyz( pyz_VelX )   

  case (IDDisturbWindSounding)  
    ! 
    !
    call MessageNotify( "M", "main", "Making Initial wind data (disturb) named Sounding..." )
    call initialdata_sounding_wind( pyz_VelX, xqz_VelY ) 
    call SetMargin_pyz( pyz_VelX )   
    call SetMargin_xqz( xqz_VelY )   

  case (IDDisturbWindZonal)  
    ! 
    !
    call MessageNotify( "M", "main", "Making Initial wind data (disturb) named Zonal..." )

    ! 値を代入
    pyz_VelX = VelX0
    xqz_VelY = VelY0
    xyr_VelZ = VelZ0

    call SetMargin_pyz( pyz_VelX )
    call SetMargin_xqz( xqz_VelY )
    call SetMargin_xyr( xyr_VelZ )

  case (IDDisturbWindTanh)  
    ! tanh 型の風速
    !
    call MessageNotify( "M", "main", "Making Initial wind data (disturb) named tanh..." )
    call initialdata_disturb_tanh(VelMean, VelMean, Zc, Zr, pyz_VelX)
    call SetMargin_pyz( pyz_VelX )   

  end select

  ! 一括して擾乱場を与える場合
  !
  select case(IDDisturb)  
  case (IDDisturbArare4)    
    call MessageNotify( "M", "main", "Making Initial data using arare4 files ...." )
    call Arare4Initfileio_Var_Get( &
      & pyz_VelX,     & ! (out)
      & xqz_VelY,     & ! (out)
      & xyr_VelZ,     & ! (out)
      & xyz_PTemp,    & ! (out)
      & xyz_Exner,    & ! (out)
      & xyzf_QMix,    & ! (out)
      & xyz_Km,       & ! (out)
      & xyz_Kh,       & ! (out)
      & xyz_CDens   )   ! (out)

  case (IDDisturbArare4mmc)    
    call MessageNotify( "M", "main", "Making Initial data using arare4 files ...." )
    call Arare4Initfileio_MMC_Var_Get( &
      & pyz_VelX,     & ! (out)
      & xqz_VelY,     & ! (out)
      & xyr_VelZ,     & ! (out)
      & xyz_PTemp,    & ! (out)
      & xyz_Exner,    & ! (out)
      & xyzf_QMix,    & ! (out)
      & xyz_Km,       & ! (out)
      & xyz_Kh,       & ! (out)
      & xyz_CDens   )   ! (out)
  end select

  !------------------------------------------------
  ! ファイル出力
  !
  call MessageNotify( "M", "main", "Output variables into netCDF file..." )

  ! リスタートファイル作成. 基本場と擾乱場を出力. 
  !
  call ReStartFileIO_Init( FlagInitData )

  call HistoryPut( 't',     Time,      rstat)
  call HistoryPut( 'VelX',  pyz_VelX,  rstat)
  call HistoryPut( 'VelY',  xqz_VelY,  rstat)
  call HistoryPut( 'VelZ',  xyr_VelZ,  rstat)
  call HistoryPut( 'Exner', xyz_Exner, rstat)
  call HistoryPut( 'PTemp', xyz_PTemp, rstat)
  call HistoryPut( 'Km',    xyz_Km,    rstat)
  call HistoryPut( 'Kh',    xyz_Kh,    rstat)
  call HistoryPut( 'QMix',  xyzf_QMix, rstat)    
  call HistoryPut( 'CDens', xyz_CDens, rstat)    

  ! 基本場のファイル出力
  !
  xyz_VPTempBZ = xyz_PTempBZ / xyz_EffMolWtBZ  ! ファイル入力のため

  call HistoryPut( 'DensBZ',     xyz_DensBZ    , rstat)
  call HistoryPut( 'ExnerBZ',    xyz_ExnerBZ   , rstat)
  call HistoryPut( 'PTempBZ',    xyz_PTempBZ   , rstat)
  call HistoryPut( 'VPTempBZ',   xyz_VPTempBZ  , rstat)
  call HistoryPut( 'VelSoundBZ', xyz_VelSoundBZ, rstat)
  call HistoryPut( 'TempBZ',     xyz_TempBZ    , rstat)
  call HistoryPut( 'PressBZ',    xyz_PressBZ   , rstat)
  call HistoryPut( 'QMixBZ',     xyzf_QMixBZ   , rstat)
  call HistoryPut( 'EffMolWtBZ', xyz_EffMolWtBZ, rstat)
  call HistoryPut( 'HumBZ',      xyzf_HumBZ    , rstat)
  
  call ReStartFileIO_Finalize

  !----------------------------------------------------
  ! MPI END
  !
  call MPIWrapperFinalize
  
contains

  subroutine ArareAlloc
    !
    !初期化として, 配列を定義し, 値としてゼロを代入する.
    !

    !暗黙の型宣言禁止
    implicit none

    !配列割り当て
    allocate( pyz_VelX(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xqz_VelY(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyr_VelZ(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_Exner(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_PTemp(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_Km(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_Kh(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_CDens(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyzf_QMix(imin:imax, jmin:jmax, kmin:kmax, ncmax) )

    allocate( xyz_DensBZ(imin:imax,jmin:jmax,kmin:kmax))
    allocate( xyz_PressBZ(imin:imax,jmin:jmax,kmin:kmax))
    allocate( xyz_ExnerBZ(imin:imax,jmin:jmax,kmin:kmax))
    allocate( xyz_TempBZ(imin:imax,jmin:jmax,kmin:kmax))
    allocate( xyz_PTempBZ(imin:imax,jmin:jmax,kmin:kmax))
    allocate( xyz_VPTempBZ(imin:imax,jmin:jmax,kmin:kmax))
    allocate( xyz_VelSoundBZ(imin:imax,jmin:jmax,kmin:kmax))
    allocate( xyzf_QMixBZ(imin:imax,jmin:jmax,kmin:kmax,ncmax))
    allocate( xyz_EffMolWtBZ(imin:imax,jmin:jmax,kmin:kmax))
    allocate( z_TempBZ(kmin:kmax))
    allocate( z_PressBZ(kmin:kmax))
    allocate( xyzf_MolFr(imin:imax, jmin:jmax, kmin:kmax, ncmax))
    allocate( zf_MolFr(kmin:kmax, ncmax))
    allocate( xyzf_QMixDivMolWt(imin:imax,jmin:jmax,kmin:kmax, ncmax))
    allocate( xyzf_HumBZ(imin:imax,jmin:jmax,kmin:kmax, ncmax))

    pyz_VelX = 0.0d0
    xqz_VelY = 0.0d0
    xyr_VelZ = 0.0d0
    xyz_PTemp = 0.0d0
    xyz_Exner = 0.0d0
    xyz_Km = 0.0d0
    xyz_Kh = 0.0d0
    xyz_CDens = 0.0d0
    xyzf_QMix = 0.0d0

    xyz_DensBZ = 0.0d0
    xyz_PressBZ = 0.0d0
    xyz_ExnerBZ = 0.0d0
    xyz_TempBZ = 0.0d0
    xyz_PTempBZ = 0.0d0
    xyz_VPTempBZ = 0.0d0
    xyz_VelSoundBZ = 0.0d0
    xyzf_QMixBZ = 0.0d0
    xyz_EffMolWtBZ  = 0.0d0
    z_TempBZ  = 0.0d0
    z_PressBZ = 0.0d0
    xyzf_MolFr = 0.0d0
    zf_MolFr = 0.0d0 
    xyzf_QMixDivMolWt = 0.0d0
    xyzf_HumBZ = 0.0d0

  end subroutine ArareAlloc
    

  subroutine DetermineBZ

    !---------------------------------------------------------------
    ! 温度・圧力
    !

    ! 3 次元配列に格納
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          xyz_TempBZ(i,j,k)  = z_TempBZ(k)  
          xyz_PressBZ(i,j,k) = z_PressBZ(k)
        end do
      end do
    end do
    
    !のり代の値を決める
    !
    call SetMargin_xyz( xyz_TempBZ )
    call SetMargin_xyz( xyz_PressBZ)

    !---------------------------------------------------------------
    ! 混合比
    !
    !水平方向には一様
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          xyzf_MolFr(i,j,k,:) = zf_MolFr(k,:)
        end do
      end do
    end do

    !のり代の値を決める
    ! 
    call SetMargin_xyzf( xyzf_MolFr )

    !気相のモル比を混合比に変換
    do s = 1, ncmax
      xyzf_QMixBZ(:,:,:,s) = xyzf_MolFr(:,:,:,s) * MolWtWet(s) / MolWtDry
    end do
    
    !  !値が小さくなりすぎないように最低値を与える
    !  where (xyzf_QMixBZ <= 1.0d-20 )
    !    xyzf_QMixBZ = 1.0d-20
    !  end where
    
    !のり代の値を決める
    !
    call SetMargin_xyzf( xyzf_QMixBZ )
    
    !---------------------------------------------------------------
    ! 分子量の効果
    !
    do s = 1, ncmax
      xyzf_QMixDivMolWt(:,:,:,s) = xyzf_QMixBZ(:,:,:,s) / MolWtWet(s)
    end do
    
    xyz_EffMolWtBZ = &
      & (1.0d0 + sum(xyzf_QMixBZ,4) ) &
      & / ( MolWtDry * ((1.0d0 / MolWtDry) + sum(xyzf_QMixDivMolWt,4)) )
    
    !のり代の値を決める  
    !
    call SetMargin_xyz( xyz_EffMolWtBZ )

    !---------------------------------------------------------------    
    ! 温位
    !
    xyz_PTempBZ = &
      & xyz_TempBZ * (PressBasis / xyz_PressBZ) ** (GasRDry / CpDry) 

    !のり代の値を決める  
    !
    call SetMargin_xyz( xyz_PTempBZ )

    
    !---------------------------------------------------------------    
    ! エクスナー関数
    !
    xyz_ExnerBZ = xyz_TempBZ / xyz_PTempBZ    
    
    !のり代の値を決める
    !
    call SetMargin_xyz( xyz_ExnerBZ )

    !---------------------------------------------------------------    
    ! 密度
    !    
    xyz_DensBZ = &
      & PressBasis * (xyz_ExnerBZ ** (CvDry / GasRDry)) &
      &  / (GasRDry * xyz_PTempBZ / xyz_EffMolWtBZ)
    
    !のり代の値を決める
    !
    call SetMargin_xyz( xyz_DensBZ )
    
    !---------------------------------------------------------------    
    ! 音速
    !
    xyz_VelSoundBZ = &
      & sqrt ( &
      &   CpDry * GasRDry * xyz_ExnerBZ * xyz_PTempBZ &
      &   / (CvDry * xyz_EffMolWtBZ) &
      & )
    
    !のり代の値を決める
    !
    call SetMargin_xyz( xyz_VelSoundBZ )

    !---------------------------------------------------------------    
    ! 湿度
    !
    if (CondNum >= 1) then 
      do s = 1, CondNum
        xyzf_HumBZ(:,:,:,IdxCG(s)) = &
          & xyz_PressBZ * xyzf_MolFr(:,:,:,IdxCG(s)) &
          & / xyz_SvapPress(SpcWetID(IdxCC(s)), xyz_TempBZ) &
          & * 100.0d0
      end do
    end if

    !----------------------------------------------------------
    ! BasicSet モジュールに値を設定
    !   
    call BasicSet_Init(                           &
      & xyz_PressBZ, xyz_ExnerBZ, xyz_TempBZ,     &
      & xyz_PTempBZ, xyz_DensBZ,  xyz_VelSoundBZ, &
      & xyzf_QMixBZ, xyz_EffMolWtBZ )
    
  end subroutine DetermineBZ


  subroutine Initialdata_init

    use dc_iounit,     only : FileOpen    

    implicit none

    integer                       :: unit     !装置番号

    character(STRING) :: FlagBasic        = ""
    character(STRING) :: FlagDisturbPTemp = ""
    character(STRING) :: FlagDisturbExner = ""
    character(STRING) :: FlagDisturbQMix  = ""
    character(STRING) :: FlagDisturbWind  = ""
    character(STRING) :: FlagDisturb      = ""

    NAMELIST /initialdata_nml/ &
      & FlagBasic,        & 
      & FlagDisturb,      &
      & FlagDisturbPTemp, &    
      & FlagDisturbExner, &    
      & FlagDisturbQMix,  &    
      & FlagDisturbWind  

    NAMELIST /initialdata_basic_nml/ &
!      & Humidity, TempTr, DHeight, HeightTr
      & Humidity, TempTr, DHeight

    NAMELIST /initialdata_disturb_gauss_nml/ &
      & PTempMax, ExnerMax, QMixMax, Xc, Xr, Yc, Yr, Zc, Zr

    NAMELIST /initialdata_disturb_cos_nml/ &
      & PTempMax, ExnerMax, QMixMax, Xc, Xr, Yc, Yr, Zc, Zr

    NAMELIST /initialdata_disturb_cone_nml/ &
      & PTempMax, ExnerMax, QMixMax, Xc, Xr, Yc, Yr, Zc, Zr

    NAMELIST /initialdata_disturb_tanh_nml/ &
      & PTempMean, PTempDel, VelMean, Zc, Zr

    NAMELIST /initialdata_disturb_random_nml/ &
!      & PTempMax, ExnerMax, QMixMax, Zpos
      & PTempMax, Zpos

    NAMELIST /initialdata_disturb_square_nml/ &
      & QMixMax, PTempMax, XposMin, YposMin, ZposMin, XposMax, YposMax, ZposMax

!    NAMELIST /initialdata_disturb_dryreg_nml/ &
!      & XposMin, YposMin, ZposMin, XposMax, YposMax, ZposMax

    NAMELIST /initialdata_disturb_zonal_nml/ &
      & VelX0, VelY0, VelZ0

    NAMELIST /initialdata_disturb_moist_nml/ &
      & Humidity

    NAMELIST /initialdata_disturb_SK1994_nml/ &
      & PTempMax, Xc, Xr, Yc, Yr

    NAMELIST /initialdata_disturb_circle_nml/ &
      & PTempMax, ExnerMax, QMixMax, Xc, Xr, Yc, Zc

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=initialdata_nml)
    close(unit)
    
    if (FlagBasic == "Arare4") then 
      IDBasic = IDBasicArare4
    elseif(FlagBasic == "Arare4mmc") then 
      IDBasic = IDBasicArare4mmc
    elseif(FlagBasic == "Yamasaki1983") then 
      IDBasic = IDBasicYamasaki1983
    elseif(FlagBasic == "Takemi2007") then 
      IDBasic = IDBasicTakemi2007
    elseif(FlagBasic == "IsoThermal") then 
      IDBasic = IDBasicIsoThermal
    elseif(FlagBasic == "Dry"       ) then 
      IDBasic = IDBasicDry
    elseif(FlagBasic == "Moist"     ) then 
      IDBasic = IDBasicMoist
    elseif(FlagBasic == "Toon2002"  ) then 
      IDBasic = IDBasicToon2002
    elseif(FlagBasic == "Sounding") then 
      IDBasic = IDBasicSounding
    elseif(FlagBasic == "Baker1998") then 
      IDBasic = IDBasicBaker1998
    elseif(FlagBasic == "SK1994") then 
      IDBasic = IDBasicSK1994
    end if

    if(FlagDisturbPTemp == "GaussXY" ) then     
      IDDisturbPTemp = IDDisturbPTempGaussXY
    elseif(FlagDisturbPTemp == "GaussXZ" ) then     
      IDDisturbPTemp = IDDisturbPTempGaussXZ
    elseif(FlagDisturbPTemp == "GaussYZ" ) then     
      IDDisturbPTemp = IDDisturbPTempGaussYZ
    elseif(FlagDisturbPTemp == "GaussXYZ") then     
      IDDisturbPTemp = IDDisturbPTempGaussXYZ
    elseif(FlagDisturbPTemp == "Random"  ) then     
      IDDisturbPTemp = IDDisturbPTempRandom
    elseif(FlagDisturbPTemp == "Square"  ) then     
      IDDisturbPTemp = IDDisturbPTempSquare
    elseif(FlagDisturbPTemp == "CosXY" ) then     
      IDDisturbPTemp = IDDisturbPTempCosXY
    elseif(FlagDisturbPTemp == "CosXZ" ) then     
      IDDisturbPTemp = IDDisturbPTempCosXZ
    elseif(FlagDisturbPTemp == "CosYZ" ) then     
      IDDisturbPTemp = IDDisturbPTempCosYZ
    elseif(FlagDisturbPTemp == "CosXYZ") then     
      IDDisturbPTemp = IDDisturbPTempCosXYZ
    elseif(FlagDisturbPTemp == "ConeXY" ) then     
      IDDisturbPTemp = IDDisturbPTempConeXY
    elseif(FlagDisturbPTemp == "ConeXZ" ) then     
      IDDisturbPTemp = IDDisturbPTempConeXZ
    elseif(FlagDisturbPTemp == "ConeYZ" ) then     
      IDDisturbPTemp = IDDisturbPTempConeYZ
    elseif(FlagDisturbPTemp == "tanh" ) then     
      IDDisturbPTemp = IDDisturbPTempTanh
    elseif(FlagDisturbPTemp == "SK1994XZ" ) then     
      IDDisturbPTemp = IDDisturbPTempSK1994XZ
    elseif(FlagDisturbPTemp == "SK1994YZ" ) then     
      IDDisturbPTemp = IDDisturbPTempSK1994YZ
    elseif(FlagDisturbPTemp == "CircleXZ" ) then     
      IDDisturbPTemp = IDDisturbPTempCircleXZ
    end if

    if(FlagDisturbExner == "GaussXY" ) then     
      IDDisturbExner = IDDisturbExnerGaussXY 
    elseif(FlagDisturbExner == "GaussXZ" ) then     
      IDDisturbExner = IDDisturbExnerGaussXZ 
    elseif(FlagDisturbExner == "GaussYZ" ) then     
      IDDisturbExner = IDDisturbExnerGaussYZ 
    elseif(FlagDisturbExner == "GaussXYZ") then    
      IDDisturbExner = IDDisturbExnerGaussXYZ  
    elseif(FlagDisturbExner == "CosXY" ) then     
      IDDisturbExner = IDDisturbExnerCosXY 
    elseif(FlagDisturbExner == "CosXZ" ) then     
      IDDisturbExner = IDDisturbExnerCosXZ 
    elseif(FlagDisturbExner == "CosYZ" ) then     
      IDDisturbExner = IDDisturbExnerCosYZ 
    elseif(FlagDisturbExner == "CosXYZ") then    
      IDDisturbExner = IDDisturbExnerCosXYZ  
    elseif(FlagDisturbExner == "ConeXY" ) then     
      IDDisturbExner = IDDisturbExnerConeXY 
    elseif(FlagDisturbExner == "ConeXZ" ) then     
      IDDisturbExner = IDDisturbExnerConeXZ 
    elseif(FlagDisturbExner == "ConeYZ" ) then     
      IDDisturbExner = IDDisturbExnerConeYZ 
    end if

    if(FlagDisturbQMix == "GaussXY"  ) then     
      IDDisturbQMix = IDDisturbQMixGaussXY
    elseif(FlagDisturbQMix == "GaussXZ"  ) then 
      IDDisturbQMix = IDDisturbQMixGaussXZ   
    elseif(FlagDisturbQMix == "GaussYZ"  ) then 
      IDDisturbQMix = IDDisturbQMixGaussYZ   
    elseif(FlagDisturbQMix == "GaussXYZ" ) then 
      IDDisturbQMix = IDDisturbQMixGaussXYZ    
    elseif(FlagDisturbQMix == "Dryreg" ) then
      IDDisturbQMix = IDDisturbQMixDryreg
    elseif(FlagDisturbQMix == "Moist"    ) then    
      IDDisturbQMix = IDDisturbQMixMoist
    elseif(FlagDisturbQMix == "CosXY"  ) then     
      IDDisturbQMix = IDDisturbQMixCosXY
    elseif(FlagDisturbQMix == "CosXZ"  ) then 
      IDDisturbQMix = IDDisturbQMixCosXZ   
    elseif(FlagDisturbQMix == "CosYZ"  ) then 
      IDDisturbQMix = IDDisturbQMixCosYZ   
    elseif(FlagDisturbQMix == "CosXYZ" ) then 
      IDDisturbQMix = IDDisturbQMixCosXYZ  
    elseif(FlagDisturbQMix == "ConeXY"  ) then     
      IDDisturbQMix = IDDisturbQMixConeXY
    elseif(FlagDisturbQMix == "ConeXZ"  ) then 
      IDDisturbQMix = IDDisturbQMixConeXZ   
    elseif(FlagDisturbQMix == "ConeYZ"  ) then 
      IDDisturbQMix = IDDisturbQMixConeYZ     
    elseif(FlagDisturbQMix == "Square") then     
      IDDisturbQMix = IDDisturbQMixSquare
    end if

    if (FlagDisturbWind == "Takemi2007") then
      IDDisturbWind = IDDisturbWindTakemi2007
    elseif (FlagDisturbWind == "Sounding") then
      IDDisturbWind = IDDisturbWindSounding
    elseif (FlagDisturbWind == "Zonal") then
      IDDisturbWind = IDDisturbWindZonal
    elseif (FlagDisturbWind == "tanh") then
      IDDisturbWind = IDDisturbWindTanh
    end if

    if(FlagDisturb == "Arare4"        ) then    
      IDDisturb = IDDisturbArare4
    elseif(FlagDisturb == "Arare4mmc" ) then    
      IDDisturb = IDDisturbArare4mmc
    end if

    !----------------------------------------------------
    ! モジュールの初期化
    !
    !   Yamasaki, baker, シンプル擾乱場は, 初期化する必要ない. 
    !
    select case (IDDisturbPTemp)
    case (IDDisturbPTempGaussXY:IDDisturbPTempGaussXYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_gauss_nml)
      close(unit)

    case (IDDisturbPTempRandom)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_random_nml)
      close(unit)

    case (IDDisturbPTempSquare)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_square_nml)
      close(unit)

    case (IDDisturbPTempCosXY:IDDisturbPTempCosXYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_cos_nml)
      close(unit)

    case (IDDisturbPTempConeXY:IDDisturbPTempConeYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_cone_nml)
      close(unit)

    case (IDDisturbPTempTanh)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_tanh_nml)
      close(unit)

    case (IDDisturbPTempSK1994XZ:IDDisturbPTempSK1994YZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_SK1994_nml)
      close(unit)

    case (IDDisturbPTempCircleXZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_circle_nml)
      close(unit)
    end select

    select case (IDDisturbExner)
    case (IDDisturbExnerGaussXY:IDDisturbExnerGaussXYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_gauss_nml)
      close(unit)

    case (IDDisturbExnerCosXY:IDDisturbExnerCosXYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_cos_nml)
      close(unit)

    case (IDDisturbExnerConeXY:IDDisturbExnerConeYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_cone_nml)
      close(unit)
    end select

    select case (IDDisturbQMix)
    case (IDDisturbQMixGaussXY:IDDisturbQMixGaussXYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_gauss_nml)
      close(unit)

    case(IDDisturbQMixDryreg)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_square_nml)
      close(unit)

    case(IDDisturbQMixMoist)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_moist_nml)
      close(unit)

    case (IDDisturbQMixCosXY:IDDisturbQMixCosXYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_cos_nml)
      close(unit)

    case (IDDisturbQMixConeXY:IDDisturbQMixConeYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_cone_nml)
      close(unit)

    case (IDDisturbQMixSquare)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_square_nml)
      close(unit)
      
    end select

    select case (IDDisturbWind)
    case (IDDisturbWindZonal)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_zonal_nml)
      close(unit)

    case (IDDisturbWindTanh)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_tanh_nml)
      close(unit)
    end select

    select case (IDBasic)
    case (IDBasicDry, IDBasicMoist)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_basic_nml)
      close(unit)
    end select

    if (FlagBasic == "Takemi2007" .OR. FlagDisturbWind == "Takemi2007") then 
      call initialdata_takemi2007_init
    end if

    if (FlagBasic == "Sounding" .OR. FlagDisturbWind == "Sounding") then 
      call initialdata_sounding_init
    end if

    if (     FlagBasic == "Arare4"        &
      & .OR. FlagDisturb == "Arare4"      &
      & .OR. FlagBasic   == "Arare4mmc"   &
      & .OR. FlagDisturb == "Arare4mmc") then 
      call Arare4InitFileio_init
    end if

    call MessageNotify( "M", "main", "FlagBasic        = %c", c1=trim(FlagBasic))
    call MessageNotify( "M", "main", "IDBasic          = %d", i=(/IDBasic/))
    call MessageNotify( "M", "main", "FlagDisturbPTemp = %c", c1=trim(FlagDisturbPTemp))
    call MessageNotify( "M", "main", "IDDisturbPTemp   = %d", i=(/IDDisturbPTemp/))
    call MessageNotify( "M", "main", "FlagDisturbExner = %c", c1=trim(FlagDisturbExner))
    call MessageNotify( "M", "main", "IDDisturbExner   = %d", i=(/IDDisturbExner/))
    call MessageNotify( "M", "main", "FlagDisturbWind  = %c", c1=trim(FlagDisturbWind))
    call MessageNotify( "M", "main", "IDDisturbWind    = %d", i=(/IDDisturbWind/))
    call MessageNotify( "M", "main", "FlagDisturbQMix  = %c", c1=trim(FlagDisturbQMix))
    call MessageNotify( "M", "main", "IDDisturbQMix    = %d", i=(/IDDisturbQMix/))
    call MessageNotify( "M", "main", "FlagDisturb      = %c", c1=trim(FlagDisturb))
    call MessageNotify( "M", "main", "IDDisturb        = %d", i=(/IDDisturb/))

  end subroutine Initialdata_init


!!!--------------------------------------------------------------
  subroutine initialdata_basic_strat(z_TempBZ, z_PressBZ)
    !
    implicit none
    
    real(DP), intent(inout) :: z_TempBZ(kmin:kmax)
    real(DP), intent(inout) :: z_PressBZ(kmin:kmax)
    real(DP)                :: HeightTr
    real(DP)                :: z_DTempDZ(kmin:kmax)
    real(DP)                :: DTempDZ
    real(DP)                :: Weight
    integer                 :: k    
    
    do k = 1, kmax
      z_DTempDZ(k) = (z_TempBZ(k) - z_TempBZ(k-1)) / dz
    end do

    ! 対流圏界面
    !
    HeightTr =  minval(z_z(1:nz), 1, z_TempBZ(1:nz) < TempTr)
    
    ! 対流圏界面より上の扱い
    !   圏界面より上は等温大気とする. 等温位大気から等温大気への遷移は 
    !   tanh を用いてなめらかにつなぐ
    do k = 2, kmax
      
      !重みつけの関数を用意. tanh を用いる
      Weight = ( tanh( (z_Z(k) - HeightTr ) / DHeight ) + 1.0d0 ) * 5.0d-1
      
      !仮値として温度を計算する. 圏界面より上では TempTr の等温大気に近づける
      z_TempBZ(k) = z_TempBZ(k) * ( 1.0d0 - Weight ) + TempTr * Weight
      
      !温度減率が断熱温度減率より小さくならないように
      DTempDZ = max( z_DTempDZ(k), (z_TempBZ(k) - z_TempBZ(k-1)) / dz )
      
      !基本場の温度を決める
      z_TempBZ(k) = z_TempBZ(k-1) + DTempDZ * dz
      
      !圧力を静水圧平衡から計算
      z_PressBZ(k) =                                          &
        &  z_PressBZ(k-1) * ( ( z_TempBZ(k-1) / z_TempBZ(k) ) &
        &    ** (Grav / ( DTempDZ * GasRDry ) ) )
    end do
    
  end subroutine initialdata_basic_strat


!!!--------------------------------------------------------------
  subroutine initialdata_basic_strat_v2(z_TempBZ, z_PressBZ)
    !
    implicit none
    
    real(DP), intent(inout) :: z_TempBZ(kmin:kmax)
    real(DP), intent(inout) :: z_PressBZ(kmin:kmax)
    real(DP)                :: z_DTempDZ(kmin:kmax)
    real(DP)                :: DTempDZ
    integer                 :: k, k1, k2

    ! 現在の温度勾配を計算. 
    !
    do k = 1, kmax
      z_DTempDZ(k) = (z_TempBZ(k) - z_TempBZ(k-1)) / dz
    end do

    ! 対流圏界面より上の扱い. 指定された温度となる高度を圏界面と見なす
    !
!    k1 = minloc( z_Z, 1, z_Z > HeightTr ) 
    k1 = maxloc( z_TempBZ, 1, z_TempBZ <= TempTr ) 

    ! 温度減率
    !
    k2 = int( DHeight / dz )
    k1 = k1 - k2 / 2          ! 指定された高度をまたぐ
    DTempDZ = z_DTempDZ(k1-1)

    do k = k1, kmax

      ! 温度源率をゼロに近づける. 
      !
      DTempDZ = min( -1.0d-14, DTempDZ - DTempDZ / k2 * ( k - k1 ) )
      
      !基本場の温度を決める
      !
      z_TempBZ(k) = z_TempBZ(k-1) + DTempDZ * dz
      
      !圧力を静水圧平衡から計算
      !
      z_PressBZ(k) =                                          &
        &  z_PressBZ(k-1) * ( ( z_TempBZ(k-1) / z_TempBZ(k) ) &
        &    ** (Grav / ( DTempDZ * GasRDry ) ) )

    end do
    
  end subroutine initialdata_basic_strat_v2


  subroutine MainInit

    implicit none

    character(STRING)  :: cfgfile = ""    
    integer, parameter :: OutputRank = 0

    ! gtool のメッセージの抑制
    ! rank0 のみ標準出力
    !
    call MessageSuppressMPI( OutputRank )

    ! MPI
    !
    call MPIWrapperInit
    
    !NAMELIST ファイル名の読み込み
    !
    call argset_init( cfgfile ) !(out)
    
    ! NAMELIST ファイル名のモジュールへの登録
    ! Loading NAMELIST file.
    !
    call NmlutilInit( cfgfile ) !(in)
    
    !格子点情報の初期化
    !  NAMELIST から情報を得て, 格子点を計算する
    !
    call gridset_init
    
    ! 化学計算ルーチンの初期化
    ! Initialization of chemical routines.
    !
    call chemcalc_init
    
    ! 軸の情報の初期化
    ! Initialization of axis variables.
    !
    call axesset_init
    
    ! 定数の情報の初期化
    ! Initialization of constant variables.
    !
    call constants_init
    
    ! 湿潤過程共有変数の初期化
    ! Initialization of common variables for moist process.
    !
    call composition_init
    
    ! I/O ファイル名の初期化
    ! Initialization of output file name. 
    !
    call fileset_init
    
    ! nml から情報を取り出す (内部サブルーチン)
    !
    call InitialData_init

    ! マージンの設定の初期化
    ! Initialization of margin
    !
    call SetMargin_Init

    !内部変数の初期化. とりあえずゼロを入れて値を確定させておく. 
    !
    call ArareAlloc  
    
  end subroutine MainInit

end program ArareInitData
