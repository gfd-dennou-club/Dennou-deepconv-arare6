!= Module DynamicsHEVI
!
! Authors::   杉山耕一朗(SUGIYAMA Ko-ichiro)
! Version::   $Id: dynamics_hevi_v3.f90,v 1.4 2014/07/08 00:55:26 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2014. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module DynamicsHEVI_v3
  !
  ! 力学コア. 
  !   タイムスプリット法を利用. 音波モードとそれ以外を別々の時間刻みで解く. 
  !   短い時間ステップの計算には HE-VI 法を利用.
  !
  ! v3: v1 と v2 を整理することで, 計算方法を切り替えられるように調整した版
  !
  ! Note: 
  !  * エクスナー関数の空間方向の離散化は 2 次精度であるため, 気圧傾度
  !    力項の計算プログラムにおいて differentiate_center4 モジュールを
  !    指定することはできないので注意.

  !モジュール読み込み
  use dc_types,   only : DP

  !暗黙の型宣言禁止
  implicit none

  !属性の指定
  private

  !変数の定義
  character(*), parameter:: module_name = 'dynamics_hevi'
                                           ! モジュールの名称.
                                           ! Module name
  real(DP), save, private :: FactorBuoyTemp    = 1.0d0
                                           !浮力 (温度の寄与) の有無
                                           !考慮しない場合は値をゼロにする.
  real(DP), save, private :: FactorBuoyMolWt   = 1.0d0
                                           !浮力 (分子量効果) の有無
                                           !考慮しない場合は値をゼロにする.
  real(DP), save, private :: FactorBuoyLoading = 1.0d0
                                           !浮力 (荷重効果) の有無
                                           !考慮しない場合は値をゼロにする.
  real(DP), save, private :: FactorDExnerDtAdv    = 0.0d0
                                           !エクスナー関数の移流の有無
                                           !考慮しない場合は値をゼロにする.
  real(DP), save, private :: FactorDExnerDtExpnd  = 0.0d0
                                           !エクスナー関数の膨張項の有無
                                           !考慮しない場合は値をゼロにする.
  real(DP), save, private :: FactorExnerFall  = 0.0d0
                                           !落下に伴う仮温位変化に伴うエクスナー関数のtendency
                                           !考慮しない場合は値をゼロにする.
  integer, save, private  :: IDAdvection               = 0
                                           ! 移流計算方法
  integer, parameter      :: IDAdvectionCenter4_std    = 1
                                           ! 1: 4 次中央差分 (微分平均モジュール利用)
  integer, parameter      :: IDAdvectionCenter4_2D_dry = 2
                                           ! 2: 4 次中央差分 (2D 版計算ルーチン)
  integer, parameter      :: IDAdvectionCenter4_2D     = 3
                                           ! 2: 4 次中央差分 (2D 版計算ルーチン)
  integer, parameter      :: IDAdvectionCenter4_3D_dry = 4
                                           ! 3: 4 次中央差分 (3D 版計算ルーチン)
  integer, parameter      :: IDAdvectionCenter4_3D      = 5
                                           ! 3: 4 次中央差分 (3D 版計算ルーチン)
  integer, save, private  :: IDAcousticmode            = 0 
                                           ! 短い時間ステップの計算方法
  integer, parameter      :: IDAcousticmode_std        = 1 
                                           ! 1: 微分平均モジュールの利用. 
  integer, parameter      :: IDAcousticmode_2D         = 2
                                           ! 2: 2D 版計算ルーチン
  integer, parameter      :: IDAcousticmode_3D         = 3
                                           ! 3: 3D 版計算ルーチン
  
  ! public 
  !
  public Dynamics_Init
  public Dynamics_Long_forcing
  public Dynamics_Short_forcing

contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subroutine Dynamics_Init
    !
    ! 力学コア 3D 版の初期化ルーチン
    !

    ! モジュール読み込み
    !
    use dc_types,   only : DP, STRING
    use dc_iounit,  only : FileOpen
    use dc_message, only : MessageNotify
    use gridset,    only : FlagCalc3D, FlagCalcMoist
    use namelist_util, &
      &             only : namelist_filename
    use advection_center4_std, &
      &             only : advection_center4_std_init
    use advection_center4_2D, &
      &             only : advection_center4_2D_init
    use advection_center4_3D, &
      &             only : advection_center4_3D_init
    use acousticmode_std, &
      &             only : acousticmode_std_init
    use acousticmode_2D, &
      &             only : acousticmode_2D_init
    use acousticmode_3D, &
      &             only : acousticmode_3D_init

    !暗黙の型宣言禁止
    !
    implicit none
    
    ! 変数の定義
    !
    real(DP)          :: AlphaSound = 5.0d-2    !音波減衰項の係数 (気象庁数値予報課報告・別冊49 より)
    real(DP)          :: AlphaNDiff = 1.0d-3    !4次の数値拡散の係数. CReSS マニュアルより
    real(DP)          :: NDiffRatio = 1.0d0     !速度に対する粘性を上げる場合は数字を 1 以上にする. 
    integer           :: unit                   !装置番号
    character(STRING) :: FlagAdvection = ""     !移流計算方法の選択
    character(STRING) :: FlagAcousticmode = ""  !音波モードの計算方法の選択    

    
    !-------------------------------------------------------------------
    ! Namelist から情報を取得する
    !
    NAMELIST /Dynamics_nml/                                 &
      & AlphaSound, AlphaNDiff, NDiffRatio,                 &
      & FactorBuoyTemp, FactorBuoyMolWt, FactorBuoyLoading, &
      & FactorDExnerDtAdv, FactorDExnerDtExpnd, FactorExnerFall,  &
      & FlagAdvection, FlagAcousticmode

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=dynamics_nml)
    close(unit)
    
    !-------------------------------------------------------------------
    ! フラグ処理
    !
    if ( FlagAdvection == "Center4" .OR. FlagAdvection == "" ) then 
      !
      ! 4 次精度中心差分. 微分平均モジュールを使わない場合. (デフォルト)
      !
      if ( FlagCalc3D ) then 
!        if ( FlagCalcTracer .OR. FlagCalcMoist ) then 
        if ( FlagCalcMoist ) then 
          IDAdvection = IDAdvectionCenter4_3D
        else
          IDAdvection = IDAdvectionCenter4_3D_dry
        end if
      else
!        if ( FlagCalcTracer .OR. FlagCalcMoist ) then 
        if ( FlagCalcMoist ) then 
          IDAdvection = IDAdvectionCenter4_2D
        else
          IDAdvection = IDAdvectionCenter4_2D_dry
        end if
      end if
    else 
      !
      ! 4 次精度中心差分 with 微分平均モジュール
      !
      IDAdvection = IDAdvectionCenter4_std
    end if
    
    if ( FlagAcousticmode == "Center2" .OR. FlagAcousticmode == "" ) then 
      !
      ! 微分平均モジュールを使わない場合. (デフォルト)
      !
      if ( FlagCalc3D ) then 
        IDAcousticmode = IDAcousticmode_3D
      else
        IDAcousticmode = IDAcousticmode_2D
      end if
    else
      !
      ! with 微分平均モジュール
      !
      IDAcousticmode = IDAcousticmode_std
    end if

    !-------------------------------------------------------------------
    ! 移流計算用の計算モジュールの初期化:
    !
    select case ( IDAdvection )

    case ( IDAdvectionCenter4_std )
      call advection_center4_std_init( AlphaNDiff, NDiffRatio )

    case ( IDAdvectionCenter4_3D_dry, IDAdvectionCenter4_3D )
      call advection_center4_3D_init( AlphaNDiff, NDiffRatio )
      
    case ( IDAdvectionCenter4_2D_dry, IDAdvectionCenter4_2D )
      call advection_center4_2D_init( AlphaNDiff, NDiffRatio )
      
    end select

    !-------------------------------------------------------------------
    ! 音波モードの計算モジュールの初期化:
    !
    select case ( IDAcousticmode )
      
    case ( IDAcousticmode_std )
      call acousticmode_std_init( AlphaSound )

    case ( IDAcousticmode_3D )
      call acousticmode_3D_init( AlphaSound )

    case ( IDAcousticmode_2D )
      call acousticmode_2D_init( AlphaSound )
      
    end select
    
    !-------------------------------------------------------------------
    ! tendency の出力
    !
    call Dynamics_Tendency_output

    !-------------------------------------------------------------------
    ! 出力
    !
    call MessageNotify( "M", module_name, "AlphaSound = %f", d=(/AlphaSound/) )
    call MessageNotify( "M", module_name, "AlphaNDiff = %f", d=(/AlphaNDiff/) )
    
    call MessageNotify( "M", module_name, "FactorBuoyTemp   = %f", d=(/FactorBuoyTemp/) )
    call MessageNotify( "M", module_name, "FactorBuoyMolWt  = %f", d=(/FactorBuoyMolWt/) )
    call MessageNotify( "M", module_name, "FactorBuoyLoading= %f", d=(/FactorBuoyLoading/) )
    call MessageNotify( "M", module_name, "FactorDExnerDtAdv   = %f", d=(/FactorDExnerDtAdv/) )
    call MessageNotify( "M", module_name, "FactorDExnerDtExpnd = %f", d=(/FactorDExnerDtExpnd/) )
    
    call MessageNotify( "M", module_name, "IDAdvection = %d", i=(/IDAdvection/) )
    call MessageNotify( "M", module_name, "IDAcousticmode = %d", i=(/IDAcousticmode/) )
 
    
  end subroutine Dynamics_Init

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Dynamics_Long_forcing(   &
    & pyz_VelXBl,  pyz_VelXNl,        & ! (in)
    & xqz_VelYBl,  xqz_VelYNl,        & ! (in)
    & xyr_VelZBl,  xyr_VelZNl,        & ! (in)
    & xyz_PTempBl, xyz_PTempNl,       & ! (in)
    & xyz_ExnerBl, xyz_ExnerNl,       & ! (in)
    & xyzf_QMixBl, xyzf_QMixNl,       & ! (in)
    & xyz_KmBl,    xyz_KmNl,          & ! (in)
    & pyz_DVelXDtNl,                  & ! (inout)
    & xqz_DVelYDtNl,                  & ! (inout)
    & xyr_DVelZDtNl,                  & ! (inout)
    & xyz_DPTempDtNl,                 & ! (inout)
    & xyz_DExnerDtNl,                 & ! (inout)
    & xyzf_DQMixDtNl,                 & ! (inout)
    & xyz_DKmDtNl                     & ! (inout)
    & )
    !
    ! 力学コア: 長い時間ステップで評価する項の計算.
    !

    ! モジュール読み込み
    !
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax, ncmax, &
                                                    !配列サイズ
      &                   nx, ny, nz,              &!物理領域の配列サイズ
      &                   FlagCalcMoist !フラグ
    use composition, &
      &             only : SpcWetSymbol
    use timeset,    only : TimeN
    use gtool_historyauto, &
      &             only : HistoryAutoPut
    use advection_center4_std, &
      &             only : advection_center4_std_main
    use advection_center4_2D, &
      &             only : advection_center4_2D_dry, advection_center4_2D_tracer
    use advection_center4_3D, &
      &             only : advection_center4_3D_dry, advection_center4_3D_tracer
    use DExnerDt,   only : xyz_DExnerDt_xyzf

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 変数の定義
    !
    real(DP), intent(in)    :: pyz_VelXBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: pyz_VelXNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZNl(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(in)    :: xyz_PTempBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_PTempNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyzf_QMixBl(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(in)    :: xyzf_QMixNl(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(in)    :: xyz_KmBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_KmNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: pyz_DVelXDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xqz_DVelYDtNl(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(inout) :: xyr_DVelZDtNl(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(inout) :: xyz_DPTempDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyzf_DQMixDtNl(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(inout) :: xyz_DKmDtNl(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)                :: pyz_DVelXDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xqz_DVelYDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_DVelZDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_DKmDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_DExnerDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_DPTempDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyzf_QMixAdv(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
    real(DP)                :: xyzf_QMixFall(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
    real(DP)                :: pyz_VelXnDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xqz_VelYnDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_VelZnDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_KmNDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_ExnerNDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_ExnerFall(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_PTempNDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyzf_QMixNDiff(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
    real(DP)                :: xyr_BuoyT(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_BuoyM(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_BuoyD(imin:imax,jmin:jmax,kmin:kmax)
    
    integer                 :: f
    
    
    !--------------------------------------------------------------------
    ! 移流項の計算
    !
    select case ( IDAdvection )
      
    case ( IDAdvectionCenter4_std ) 

      call advection_center4_std_main(     &
        & pyz_VelXBl,   pyz_VelXNl,        & ! (in)
        & xqz_VelYBl,   xqz_VelYNl,        & ! (in)
        & xyr_VelZBl,   xyr_VelZNl,        & ! (in)
        & xyz_PTempBl,  xyz_PTempNl,       & ! (in)
        & xyz_ExnerBl,  xyz_ExnerNl,       & ! (in)
        & xyzf_QMixBl,  xyzf_QMixNl,       & ! (in)
        & xyz_KmBl,     xyz_KmNl,          & ! (in)
        & pyz_DVelXDtAdv,  pyz_VelXnDiff,     & !(out)
        & xqz_DVelYDtAdv,  xqz_VelYnDiff,     & !(out)
        & xyr_DVelZDtAdv,  xyr_VelZnDiff,     & !(out)
        & xyz_DPTempDtAdv, xyz_PTempNDiff,    & !(out)
        & xyz_DExnerDtAdv, xyz_ExnerNDiff,    & !(out)
        & xyzf_QMixAdv, xyzf_QMixNDiff,    & !(out)
        & xyz_DKmDtAdv,    xyz_KmNDiff        & !(out)
        & )
      
    case ( IDAdvectionCenter4_3D_dry ) 
      
      call advection_center4_3D_dry(       &
        & pyz_VelXBl,   pyz_VelXNl,        & ! (in)
        & xqz_VelYBl,   xqz_VelYNl,        & ! (in)
        & xyr_VelZBl,   xyr_VelZNl,        & ! (in)
        & xyz_PTempBl,  xyz_PTempNl,       & ! (in)
        & xyz_ExnerBl,  xyz_ExnerNl,       & ! (in)
        & xyz_KmBl,     xyz_KmNl,          & ! (in)
        & pyz_DVelXDtAdv,  pyz_VelXnDiff,     & !(out)
        & xqz_DVelYDtAdv,  xqz_VelYnDiff,     & !(out)
        & xyr_DVelZDtAdv,  xyr_VelZnDiff,     & !(out)
        & xyz_DPTempDtAdv, xyz_PTempNDiff,    & !(out)
        & xyz_DExnerDtAdv, xyz_ExnerNDiff,    & !(out)
        & xyz_DKmDtAdv,    xyz_KmNDiff        & !(out)
        & )

    case ( IDAdvectionCenter4_3D ) 
      
      call advection_center4_3D_dry(       &
        & pyz_VelXBl,   pyz_VelXNl,        & ! (in)
        & xqz_VelYBl,   xqz_VelYNl,        & ! (in)
        & xyr_VelZBl,   xyr_VelZNl,        & ! (in)
        & xyz_PTempBl,  xyz_PTempNl,       & ! (in)
        & xyz_ExnerBl,  xyz_ExnerNl,       & ! (in)
        & xyz_KmBl,     xyz_KmNl,          & ! (in)
        & pyz_DVelXDtAdv,  pyz_VelXnDiff,     & !(out)
        & xqz_DVelYDtAdv,  xqz_VelYnDiff,     & !(out)
        & xyr_DVelZDtAdv,  xyr_VelZnDiff,     & !(out)
        & xyz_DPTempDtAdv, xyz_PTempNDiff,    & !(out)
        & xyz_DExnerDtAdv, xyz_ExnerNDiff,    & !(out)
        & xyz_DKmDtAdv,    xyz_KmNDiff        & !(out)
        & )

      call advection_center4_3D_tracer(       &
        & pyz_VelXNl, xqz_VelYNl, xyr_VelZNl, & ! (in)
        & xyzf_QMixBl,  xyzf_QMixNl,          & ! (in)
        & xyzf_QMixAdv, xyzf_QMixNDiff        & !(out)
        & )
      
    case ( IDAdvectionCenter4_2D_dry ) 
      
      ! Y 方向は tendency は零. 
      !
      xqz_DVelYDtAdv   = 0.0d0
      xqz_VelYnDiff = 0.0d0

      call advection_center4_2D_dry(       &
        & pyz_VelXBl,   pyz_VelXNl,        & ! (in)
        & xyr_VelZBl,   xyr_VelZNl,        & ! (in)
        & xyz_PTempBl,  xyz_PTempNl,       & ! (in)
        & xyz_ExnerBl,  xyz_ExnerNl,       & ! (in)
        & xyz_KmBl,     xyz_KmNl,          & ! (in)
        & pyz_DVelXDtAdv,  pyz_VelXnDiff,     & !(out)
        & xyr_DVelZDtAdv,  xyr_VelZnDiff,     & !(out)
        & xyz_DPTempDtAdv, xyz_PTempNDiff,    & !(out)
        & xyz_DExnerDtAdv, xyz_ExnerNDiff,    & !(out)
        & xyz_DKmDtAdv,    xyz_KmNDiff        & !(out)
        & )

    case ( IDAdvectionCenter4_2D ) 
      
      ! Y 方向は tendency は零. 
      !
      xqz_DVelYDtAdv   = 0.0d0
      xqz_VelYnDiff = 0.0d0
      
      call advection_center4_2D_dry(       &
        & pyz_VelXBl,   pyz_VelXNl,        & ! (in)
        & xyr_VelZBl,   xyr_VelZNl,        & ! (in)
        & xyz_PTempBl,  xyz_PTempNl,       & ! (in)
        & xyz_ExnerBl,  xyz_ExnerNl,       & ! (in)
        & xyz_KmBl,     xyz_KmNl,          & ! (in)
        & pyz_DVelXDtAdv,  pyz_VelXnDiff,     & !(out)
        & xyr_DVelZDtAdv,  xyr_VelZnDiff,     & !(out)
        & xyz_DPTempDtAdv, xyz_PTempNDiff,    & !(out)
        & xyz_DExnerDtAdv, xyz_ExnerNDiff,    & !(out)
        & xyz_DKmDtAdv,    xyz_KmNDiff        & !(out)
        & )

      call advection_center4_2D_tracer(    &
        & pyz_VelXNl, xyr_VelZNl,          & ! (in)
        & xyzf_QMixBl,  xyzf_QMixNl,       & ! (in)
        & xyzf_QMixAdv, xyzf_QMixNDiff     & !(out)
        & )

    end select

    !--------------------------------------------------------------------
    ! tendency の更新
    !

    ! 拡散係数 
    !
    xyz_DKmDtNl = xyz_DKmDtNl + ( xyz_KmNDiff + xyz_DKmDtAdv )
    
    ! output
    !
    call HistoryAutoPut(TimeN, 'DKmDtAdv',  xyz_DKmDtAdv(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DKmDtDiff', xyz_KmNDiff(1:nx,1:ny,1:nz) )
       
    ! 温位
    !
    xyz_DPTempDtNl = xyz_DPTempDtNl + ( xyz_PTempNDiff + xyz_DPTempDtAdv )
    
    ! output
    !
    call HistoryAutoPut(TimeN, 'DPTempDtAdv',  xyz_DPTempDtAdv(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DPTempDtDiff', xyz_PTempNDiff(1:nx,1:ny,1:nz) )

    ! 混合比
    !
    call QMixfall  ! 落下項

    xyzf_DQMixDtNl = xyzf_DQMixDtNl + ( xyzf_QMixNDiff + xyzf_QMixAdv + xyzf_QMixFall )
    
    ! output
    !
    do f = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtAdv',   &
        & xyzf_QMixAdv(1:nx,1:ny,1:nz,f))
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtDiff', &
        & xyzf_QMixNDiff(1:nx,1:ny,1:nz,f))
    end do

    ! エクスナー関数
    !
    xyz_ExnerFall   = xyz_DExnerDt_xyzf( xyzf_QMixFall ) * FactorExnerFall
    xyz_ExnerNDiff  = xyz_ExnerNDiff   * FactorDExnerDtAdv
    xyz_DExnerDtAdv = xyz_DExnerDtAdv  * FactorDExnerDtAdv
    xyz_DExnerDtNl  = xyz_DExnerDtNl + ( xyz_ExnerFall+ xyz_ExnerNDiff + xyz_DExnerDtAdv )
   
    ! output
    !
    call HistoryAutoPut(TimeN, 'DExnerDtAdv',  xyz_DExnerDtAdv(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DExnerDtDiff', xyz_ExnerNDiff(1:nx,1:ny,1:nz) )
    
    ! VelX 
    !
    pyz_DVelXDtNl = pyz_DVelXDtNl + ( pyz_VelXNDiff + pyz_DVelXDtAdv )
    
    call HistoryAutoPut(TimeN, 'DVelXDtAdv',  pyz_DVelXDtAdv(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DVelXDtDiff', pyz_VelXnDiff(1:nx,1:ny,1:nz) )
    
    ! VelY
    !
    xqz_DVelYDtNl = xqz_DVelYDtNl + ( xqz_VelYNDiff + xqz_DVelYDtAdv )
    
    call HistoryAutoPut(TimeN, 'DVelYDtAdv',  xqz_DVelYDtAdv(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DVelYDtDiff', xqz_VelYnDiff(1:nx,1:ny,1:nz) )
    
    ! VelZ
    !
    call BuoyancyLong_xyr       ! 浮力項
    
    xyr_DVelZDtNl =                            &
      &   xyr_DVelZDtNl                        &
      & + (                                    &
      &      xyr_VelZnDiff                     & 
      &    + xyr_DVelZDtAdv                       &
      &    + (                                 &
      &        + xyr_BuoyM * FactorBuoyMolWt   &
      &        + xyr_BuoyD * FactorBuoyLoading &
      &        + xyr_BuoyT * FactorBuoyTemp    &
      &      )                                 &
      &   )
    
    call HistoryAutoPut(TimeN, 'DVelZDtAdv',   xyr_DVelZDtAdv(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DVelZDtDiff',  xyr_VelZnDiff(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DVelZDtBuoyT', xyr_BuoyT(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DVelZDtBuoyM', xyr_BuoyM(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DVelZDtBuoyD', xyr_BuoyD(1:nx,1:ny,1:nz) )
    
  contains
    
    subroutine BuoyancyLong_xyr
      !
      ! 浮力項の計算
      !
    
      use composition,  only: GasNum,       &! 
        &                     IdxG,         &!
        &                     MolWtWet       ! 湿潤成分の分子量
      use constants,   only : MolWtDry,     &! 乾燥成分の分子量
        &                     Grav           ! 重力加速度
      use basicset,    only : xyz_PTempBZ,  &! 基本場の温位
        &                     xyr_QMixBZ,   &! 基本場の混合比
        &                     xyr_QMixBZPerMolWt

      ! 暗黙の型宣言禁止
      !
      implicit none
      
      ! 作業変数
      !
      real(DP)    :: xyzf_QMixPerMolWt(imin:imax,jmin:jmax,kmin:kmax, 1:GasNum)
      real(DP)    :: xyz_QMixPerMolWtSum(imin:imax,jmin:jmax,kmin:kmax)
      real(DP)    :: xyz_QMixNlSumGas(imin:imax,jmin:jmax,kmin:kmax)
      real(DP)    :: xyz_QMixNlSumCnd(imin:imax,jmin:jmax,kmin:kmax)
      real(DP)    :: xyz_BuoyT(imin:imax,jmin:jmax,kmin:kmax)
      integer     :: i, j, k, f, n
      
      ! Buoyancy due to temperature disturbunce
      !
      xyz_BuoyT = Grav * xyz_PTempNl / xyz_PTempBZ
      
      ! 格子位置の変換
      !
      do k = kmin, kmax-1
        xyr_BuoyT(:,:,k) =                           &
          & (                                        &
          &   xyz_BuoyT(:,:,k+1) + xyz_BuoyT(:,:,k)  &
          & ) * 5.0d-1 
      end do

      !穴埋め
      xyr_BuoyT(:,:,kmax) = 0.0d0

      ! 乾燥の場合は BuoyD, BuoyM は零. 
      !
      if ( .NOT. FlagCalcMoist ) then 

        xyr_BuoyM = 0.0d0
        xyr_BuoyD = 00d0

        ! サブルーチンを抜ける
        ! 
        return
      end if
           
      ! Buoyancy due to molecular weight
      !
      do f = 1, GasNum
        n = IdxG(f)
        xyzf_QMixPerMolWt(:,:,:,f) = xyzf_QMixNl(:,:,:,n) / MolWtWet(n)
      end do

      xyz_QMixPerMolWtSum = sum(xyzf_QMixPerMolWt, 4) 
      xyz_QMixNlSumGas    = sum(xyzf_QMixNl(:,:,:,1:GasNum), 4)
      
      do k = kmin, kmax-1
        do j = jmin, jmax
          do i = imin, imax

            xyr_BuoyM(i,j,k) =                                         &
              & + Grav                                                 &
              &   * (                                                  &
              &         xyz_QMixPerMolWtSum(i,j,k+1)                   &
              &       + xyz_QMixPerMolWtSum(i,j,k)                     &
              &     ) * 5.0d-1                                         &
              &   / ( 1.0d0 / MolWtDry + xyr_QMixBZPerMolWt(i,j,k) )   &
              & - Grav                                                 &
              &   * (                                                  &
              &         xyz_QMixNlSumGas(i,j,k+1)                      &
              &       + xyz_QMixNlSumGas(i,j,k)                        &
              &     ) * 5.0d-1                                         &
              &   / ( 1.0d0 + xyr_QMixBZ(i,j,k) ) 
            
          end do
        end do
      end do

      ! 穴埋め
      !
      xyr_BuoyM(:,:,kmax) = 0.0d0
      
      ! Buoyancy due to loading
      !
      xyz_QMixNlSumCnd = sum( xyzf_QMixNl(:,:,:,GasNum+1:ncmax), 4 ) 
      
      do k = kmin, kmax-1
        do j = jmin, jmax
          do i = imin, imax
            
            xyr_BuoyD(i,j,k) =                       &
              & - Grav                               &
              &   * (                                &
              &         xyz_QMixNlSumCnd(i,j,k+1)    &
              &       + xyz_QMixNlSumCnd(i,j,k)      &
              &     ) * 5.0d-1                       &
              &   / ( 1.0d0 + xyr_QMixBZ(i,j,k) )
            
          end do
        end do
      end do

      ! 穴埋め
      !
      xyr_BuoyD(:,:,kmax) = 0.0d0
      
    end subroutine BuoyancyLong_xyr


    subroutine QmixFall
      !
      ! 雨粒の落下による移流を求める. 
      ! 

      ! モジュール呼び出し
      !
      use average,          only : xyr_xyz
      use differentiate_center4, &
        &                   only : xyz_dz_xyr
      use constants,        only : FactorJ
      use composition,      only : IdxR, RainNum
      use basicset,         only : xyzf_QMixBZ,  &!基本場の混合比
        &                          xyz_DensBZ     !基本場の密度
      
      !暗黙の型宣言禁止
      !
      implicit none
      
      !変数定義
      !
      real(DP)  :: xyzf_QMixAll(imin:imax,jmin:jmax,kmin:kmax,ncmax)
                                                 !蒸気混合比(擾乱 + 平均場)
      real(DP)  :: xyrf_QMixFallFlux(imin:imax,jmin:jmax,kmin:kmax,ncmax)
                                                 !雨粒の落下効果
      real(DP)  :: xyz_VelZFall(imin:imax,jmin:jmax,kmin:kmax)
                                                 !雨粒落下速度
      integer   :: s, iR

      ! 初期化
      !
      xyzf_QMixAll = max( 0.0d0, xyzf_QMixBl + xyzf_QMixBZ )
      xyrf_QMixFallFlux = 0.0d0
      xyzf_QMixFall = 0.0d0
      xyz_VelZFall = 0.0d0

      ! 落下項の計算
      !
      do s = 1, RainNum

        iR = IdxR(s)

        !雨粒終端速度
        xyz_VelZFall = - 12.2d0 * FactorJ * ( xyzf_QMixAll(:,:,:,iR) ** 0.125d0 )
        
        ! フラックスの計算
        !
        xyrf_QMixFallFlux(:,:,:,iR) =                                &
          &  xyr_xyz (                                               &
          &    xyz_DensBZ * xyzf_QMixAll(:,:,:,iR) * xyz_VelZFall    &
          &  )

        ! 上端のフラックスはゼロ
        !        
        xyrf_QMixFallFlux(:,:,nz,iR) = 0.0d0

        ! 雨粒落下による時間変化 
        !        
        xyzf_QMixFall(:,:,:,iR) =                                      &
          &  - xyz_dz_xyr( xyrf_QMixFallFlux(:,:,:,iR) ) / xyz_DensBZ  
        
        call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(iR))//'DtFall', &
          & xyzf_QMixFall(1:nx,1:ny,1:nz,iR))
        
      end do
            
    end subroutine QmixFall

  end subroutine Dynamics_Long_forcing
  
  
!!!---------------------------------------------------------------!!!
  
  subroutine Dynamics_Short_forcing(  &
        &  pyz_VelXNs,          & ! (in)
        &  xqz_VelYNs,          & ! (in)
        &  xyr_VelZNs,          & ! (in)
        &  xyz_ExnerNs,         & ! (in)
        &  pyz_DVelXDtNl,       & ! (in)
        &  xqz_DVelYDtNl,       & ! (in)
        &  xyr_DVelZDtNl,       & ! (in)
        &  xyz_DExnerDtNl,      & ! (in)
        &  xyz_DExnerDtNs,      & ! (in)
        &  pyz_VelXAs,          & ! (out)
        &  xqz_VelYAs,          & ! (out)
        &  xyr_VelZAs,          & ! (out)
        &  xyz_ExnerAs          & ! (out)
        & )
    !
    ! 短い時間ステップで評価する項の計算.
    !
    

    ! モジュール読み込み
    !
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax, &
                                             !配列サイズ
      &                  nx, ny, nz          !物理領域の大きさ
    use timeset,  only : DelTimeShort, TimeN
    use setmargin,only : SetMargin_xyzf, SetMargin_xyz, &
      &                  SetMargin_pyz, SetMargin_xqz, SetMargin_xyr
    use gtool_historyauto, &
      &           only : HistoryAutoPut
    use acousticmode_std, &
      &           only : acousticmode_std_exp, acousticmode_std_imp
    use acousticmode_2D, &
      &           only : acousticmode_2D_exp, acousticmode_2D_imp
    use acousticmode_3D, &
      &           only : acousticmode_3D_exp, acousticmode_3D_imp
    use constants,only : CvDry,    &! 乾燥成分の比熱
      &                  GasRDry
    
    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 変数の定義
    !
    real(DP), intent(in)    :: pyz_VelXNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: pyz_DVelXDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_DVelYDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_DVelZDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_DExnerDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDtNs(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(out)   :: pyz_VelXAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xqz_VelYAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyr_VelZAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_ExnerAs(imin:imax,jmin:jmax,kmin:kmax)
    
    real(DP) :: pyz_DVelXDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xqz_DVelYDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_DVelZDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_VelDivNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_Expnd(imin:imax,jmin:jmax,kmin:kmax)

    ! subroutine の引数となる配列は, 物理領域の大きさに.
    !
    real(DP) :: pyz_PGrad(1:nx,1:ny,1:nz)
    real(DP) :: xqz_PGrad(1:nx,1:ny,1:nz)
    real(DP) :: xyr_PGrad(1:nx,1:ny,1:nz)
    real(DP) :: pyz_SWF(1:nx,1:ny,1:nz)
    real(DP) :: xqz_SWF(1:nx,1:ny,1:nz)
    real(DP) :: xyr_SWF(1:nx,1:ny,1:nz)

    !------------------------------------------------------------
    ! 初期化
    !
    pyz_DVelXDtNs = 0.0d0
    xqz_DVelYDtNs = 0.0d0
    xyr_DVelZDtNs = 0.0d0
    
    !------------------------------------------------------------
    ! 水平方向: 陽解法
    !
    select case ( IDAcousticmode )
      
    case ( IDAcousticmode_std )
      call acousticmode_std_exp(              &
        & pyz_VelXNs, xqz_VelYNs, xyr_VelZNs, & !(IN)
        & xyz_ExnerNs,                        & !(IN)
        & xyz_VelDivNs,                       & !(OUT)
        & pyz_PGrad, pyz_SWF,                 & !(OUT)
        & xqz_PGrad, xqz_SWF                  & !(OUT)
        & )
      
    case ( IDAcousticmode_3D )
      call acousticmode_3d_exp(               &
        & pyz_VelXNs, xqz_VelYNs, xyr_VelZNs, & !(IN)
        & xyz_ExnerNs,                        & !(IN)
        & xyz_VelDivNs,                       & !(IN)
        & pyz_PGrad, pyz_SWF,                 & !(OUT)
        & xqz_PGrad, xqz_SWF                  & !(OUT)
        & )
      
    case ( IDAcousticmode_2D )
      xqz_PGrad = 0.0d0
      xqz_SWF   = 0.0d0

      call acousticmode_2d_exp(               &
        & pyz_VelXNs, xyr_VelZNs,             & !(IN)
        & xyz_ExnerNs,                        & !(IN)
        & xyz_VelDivNs,                       & !(IN)
        & pyz_PGrad, pyz_SWF                  & !(OUT)
        & )
      
    end select

    ! tendency 
    !
    pyz_DVelXDtNs(1:nx,1:ny,1:nz) = pyz_PGrad + pyz_SWF
    xqz_DVelYDtNs(1:nx,1:ny,1:nz) = xqz_PGrad + xqz_SWF
    
    ! 値の保管
    !
!    call HistoryAutoPut(TimeN, 'VelDiv', xyz_VelDivNs(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'VelDiv', pyz_VelXNs(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'VelDiv', xyr_VelZNs(1:nx,1:ny,1:nz))

    call HistoryAutoPut(TimeN, 'DVelXDtPGrad', pyz_PGrad )
    call HistoryAutoPut(TimeN, 'DVelXDtSWF',   pyz_SWF   )
    call HistoryAutoPut(TimeN, 'DVelYDtPGrad', xqz_PGrad )
    call HistoryAutoPut(TimeN, 'DVelYDtSWF',   xqz_SWF   )

    ! Time integration
    !
    pyz_VelXAs = pyz_VelXNs + DelTimeShort * ( pyz_DVelXDtNl + pyz_DVelXDtNs )
    xqz_VelYAs = xqz_VelYNs + DelTimeShort * ( xqz_DVelYDtNl + xqz_DVelYDtNs )

    ! Set Margin
    !
    call SetMargin_pyz( pyz_VelXAs ) ! (inout)
    call SetMargin_xqz( xqz_VelYAs ) ! (inout)
    

    !------------------------------------------------------------
    ! 鉛直方向: 陰解法
    !  凝結項を短い時間ステップで評価する場合を考えて, 長い時間ステップと
    !  短い時間ステップで評価された tendency を足し合わせたものを引数に. 
    !

    ! 短い時間ステップで評価する圧力の式の tendency
    !
    xyz_Expnd = GasRDry * xyz_ExnerNs * xyz_VelDivNs / CvDry   &
      &         * FactorDExnerDtExpnd

    ! tendency の合計
    !
    xyz_DExnerDtNs = xyz_DExnerDtNs + xyz_Expnd

    ! 値の保管
    !
    call HistoryAutoPut(TimeN, 'DExnerDtExpnd', xyz_Expnd(1:nx,1:ny,1:nz))   

    select case ( IDAcousticmode )
      
    case ( IDAcousticmode_std )
      
      call acousticmode_std_imp(                             &
        & pyz_VelXAs, xqz_VelYAs, xyr_VelZNs, xyz_VelDivNs,  & !(IN)
        & xyz_ExnerNs,                                       & !(IN)
        & xyr_DVelZDtNl, xyz_DExnerDtNl, xyz_DExnerDtNs,     & !(IN)
        & xyz_ExnerAs,                                       & !(OUT)
        & xyr_PGrad, xyr_SWF                                 & !(OUT)
        & )
      
    case ( IDAcousticmode_3D )
      
      call acousticmode_3D_imp(                              &
        & pyz_VelXAs, xqz_VelYAs, xyr_VelZNs, xyz_VelDivNs,  & !(IN)
        & xyz_ExnerNs,                                       & !(IN)
        & xyr_DVelZDtNl, xyz_DExnerDtNl, xyz_DExnerDtNs,     & !(IN)
        & xyz_ExnerAs,                                       & !(OUT)
        & xyr_PGrad, xyr_SWF                                 & !(OUT)
        & )
      
    case ( IDAcousticmode_2D )

      call acousticmode_2D_imp(                              &
        & pyz_VelXAs, xyr_VelZNs, xyz_VelDivNs,              & !(IN)
        & xyz_ExnerNs,                                       & !(IN)
        & xyr_DVelZDtNl, xyz_DExnerDtNl, xyz_DExnerDtNs,     & !(IN)
        & xyz_ExnerAs,                                       & !(OUT)
        & xyr_PGrad, xyr_SWF                                 & !(OUT)
        & )
      
    end select
    
    ! tendency
    !
    xyr_DVelZDtNs(1:nx,1:ny,1:nz) = xyr_PGrad + xyr_SWF
    
    ! Time integration
    !
    xyr_VelZAs = xyr_VelZNs + DelTimeShort * ( xyr_DVelZDtNl + xyr_DVelZDtNs )

    ! Set Margin
    !
    call SetMargin_xyz( xyz_ExnerAs ) ! (inout)
    call SetMargin_xyr( xyr_VelZAs )  ! (inout)

    ! 値の保管
    !
    call HistoryAutoPut(TimeN, 'DVelZDtPGrad',  xyr_PGrad )
    call HistoryAutoPut(TimeN, 'DVelZDtSWF',    xyr_SWF   )
    
  end subroutine Dynamics_Short_forcing
  
!!!------------------------------------------------------------------------!!!

  subroutine Dynamics_Tendency_Output
    !
    ! 出力の設定
    ! 

    use gtool_historyauto, only : HistoryAutoAddVariable
    use gridset,           only : ncmax             ! 物質数
    use composition,       only : SpcWetSymbol, &
      &                           IdxR, RainNum
    
    implicit none
   
    integer :: f, iR

    call HistoryAutoAddVariable(                              &
      & varname='DPTempDtAdv',                                   &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Advection term of potential temperature',   &
      & units='K.s-1',                                        &
      & xtype='double')
    
    call HistoryAutoAddVariable(                                      &
      & varname='DPTempDtDiff',                                          &
      & dims=(/'x','y','z','t'/),                                     &
      & longname='Numerical diffusion term of potential temperature', &
      & units='K.s-1',                                                &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DExnerDtAdv',                                   &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Advection term of exner function',          &
      & units='s-1',                                          &
      & xtype='double')
    
    call HistoryAutoAddVariable(                               &
      & varname='DExnerDtDiff',                                   &
      & dims=(/'x','y','z','t'/),                              &
      & longname='Numerical diffusion term of exner function', &
      & units='s-1',                                           & 
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DExnerDtExpnd',                                 &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Expanding term of exner function',          &
      & units='s-1',                                          &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='CDensAdv',                                   &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Advection term of cloud density',           &
      & units='kg.m-3.s-1',                                   & 
      & xtype='double')
    
    call HistoryAutoAddVariable(                              &
      & varname='CDensDiff',                                  &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Numerical diffusion term of cloud density', &
      & units='kg.m-3.s-1',                                   &
      & xtype='double')

    do f = 1, ncmax
      call HistoryAutoAddVariable(                            &
        & varname='D'//trim(SpcWetSymbol(f))//'DtAdv',              &
        & dims=(/'x','y','z','t'/),                           &
        & longname='Advection term of '                       &
        &           //trim(SpcWetSymbol(f))//' mixing ratio', &
        & units='kg.kg-1.s-1',                                &
        & xtype='double')
      
      call HistoryAutoAddVariable(                            &
        & varname='D'//trim(SpcWetSymbol(f))//'DtDiff',             &
        & dims=(/'x','y','z','t'/),                           &
        & longname='Diffusion term of '                       &
        &           //trim(SpcWetSymbol(f))//' mixing ratio', &
        & units='kg.kg-1.s-1',                                &
        & xtype='double')
    end do

    do f = 1, RainNum
      iR = IdxR(f)
      call HistoryAutoAddVariable(                             &
        & varname='D'//trim(SpcWetSymbol(iR))//'DtFall',             &
        & dims=(/'x','y','z','t'/),                            &
        & longname='Fall term of '                             &
        &           //trim(SpcWetSymbol(iR))//' mixing ratio', &
        & units='kg.kg-1.s-1',                                 &
        & xtype='double')

    end do

    call HistoryAutoAddVariable(                              &
      & varname='DVelXDtAdv',                                    &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Advection term of velocity (x)',            &
      & units='m.s-2',                                        &
      & xtype='double')
    
    call HistoryAutoAddVariable(                              &
      & varname='DVelXDtDiff',                                   &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Numerical diffusion term of velocity (x)',  &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelXDtPGrad',                                  &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Pressure gradient term of velocity (x)',    &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelXDtSWF',                                    &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Filter for acoustic mode (x)',              &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelYDtAdv',                                    &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Advection term of velocity (y)',            &
      & units='m.s-2',                                        &
      & xtype='double')
    
    call HistoryAutoAddVariable(                              &
      & varname='DVelYDtDiff',                                   &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Numerical diffusion term of velocity (y)',  &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelYDtPGrad',                                  &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Pressure gradient term of velocity (y)',    &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelYDtSWF',                                    &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Filter for acoustic mode (y)',              &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelZDtAdv',                                    &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Advection term of velocity (z)',            &
      & units='m.s-2',                                        &
      & xtype='double')
    
    call HistoryAutoAddVariable(                              &
      & varname='DVelZDtDiff',                                   &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Numerical diffusion term of Velocity (z)',  &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelZDtBuoyT',                                  &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Buoyancy (Temperature)',                    &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelZDtBuoyM',                                  &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Buoyancy (MolWt)',                          &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelZDtBuoyD',                                  &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Buoyancy (Drag)',                           &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelZDtPGrad',                                  &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Pressure gradient term of velocity (z)',    &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelZDtSWF',                                    &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Filter for acoustic mode (z)',              &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='VelDiv',                                     &
      & dims=(/'x','y','z','t'/),                             &
      & longname='velocity divergence',                       &
      & units='s-2',                                          &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DKmDtAdv',                                      &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Advection of Km',                           &
      & units='s-1',                                          &
      & xtype='double')
    
    call HistoryAutoAddVariable(                              &
      & varname='DKmDtDiff',                                     &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Diffusion term of Km',                      &
      & units='s-1',                                          &
      & xtype='double')
    
!    call HistoryAutoAddVariable(  &
!      & varname='D', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='matrix D',  &
!      & units='s-1',    &
!      & xtype='double')

!    call HistoryAutoAddVariable(  &
!      & varname='E', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='matrix E',  &
!      & units='s-1',    &
!      & xtype='double')

!    call HistoryAutoAddVariable(  &
!      & varname='F', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='matrix F',  &
!      & units='s-1',    &
!      & xtype='double')

!    call HistoryAutoAddVariable(  &
!      & varname='ENs', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='DExnerDtNs',  &
!      & units='s-1',    &
!      & xtype='double')

!    call HistoryAutoAddVariable(  &
!      & varname='ENl', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='DExnerDtNl',  &
!      & units='s-1',    &
!      & xtype='double')

  end subroutine Dynamics_Tendency_Output
  
  
end module DynamicsHEVI_v3

