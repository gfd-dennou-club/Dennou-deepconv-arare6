!= Module Turbulence_kw1978
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA, Masatsugu
! Version::   $Id: turbulence_kw1978.f90,v 1.26 2014/06/07 17:34:27 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Turbulence_kw1978_v1
  !
  ! Klemp and Wilhelmson (1978) の乱流過程
  ! 乱流エネルギーの時間発展方程式を解くことで乱流拡散係数を決める
  !

  !モジュール読み込み 
  use dc_types, only: DP, STRING

  !暗黙の型宣言禁止
  implicit none

  !関数を public に設定
  public Turbulence_KW1978_Init
  public Turbulence_KW1978_Forcing

  !変数定義
  real(DP), save, private :: Cm     = 2.0d-1          !乱流エネルギー診断式の係数 
  real(DP), save, private :: MixLen = 0.0d0           !平均混合距離
  real(DP), save, public  :: KmMax  = 0.0d0           !乱流拡散係数の最大値

  real(DP), save, private :: FactorDExnerDtTurb = 1.0d0 !乱流拡散項を考慮するかのスイッチ
  logical,  save, private :: FlagArare4  = .true.       !基本場の密度を考慮するかのスイッチ

  real(DP), allocatable, save,private :: xyr_DPTempBZDz(:,:,:) 
                                                  !基本場の鉛直微分
  real(DP), allocatable, save,private :: xyrf_DQMixBZDz(:,:,:,:)
                                                  !基本場の鉛直微分

  character(*), parameter:: module_name = 'turbulence_kw1978'

contains

!!!------------------------------------------------------------------------!!!
  subroutine turbulence_kw1978_init
    !
    ! Turbulence モジュールの初期化ルーチン
    ! 

    !モジュール読み込み 
    !    
    use gridset,       only : imin, imax, jmin, jmax, kmin, kmax, ncmax
    use axesset,       only : dx,            &! x 方向の格子点間隔
      &                       dy,            &! y 方向の格子点間隔
      &                       dz              ! z 方向の格子点間隔
    use dc_iounit,     only : FileOpen
    use dc_message,    only : MessageNotify
    use timeset,       only : DelTimeLong
    use gridset,       only : FlagCalc3D 
    use namelist_util, only : namelist_filename
    use differentiate_center2, &
      &                only : xyr_dz_xyz
    use basicset,      only : xyz_PTempBZ,   &! 基本場の温位
      &                       xyzf_QMixBZ

    !暗黙の型宣言禁止
    !
    implicit none

    !作業変数
    !
    integer :: unit, f

    !-------------------------------------------------------------------
    ! NAMELIST から情報を取得
    !
    NAMELIST /turbulence_kw1978_nml/     &
      & Cm, KmMax,                       &
      & FactorDExnerDtTurb, FlagArare4

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=turbulence_kw1978_nml)
    close(unit)

    !-------------------------------------------------------------------
    ! 混合距離
    ! 2 次元計算の場合には DelY に依存しないようにするために if 文を利用.
    ! 
    if ( FlagCalc3D ) then 
      MixLen = ( dx * dy * dz ) ** (1.0d0 / 3.0d0)
    else
      MixLen = sqrt( dx * dz )
    end if
    
    !-------------------------------------------------------------------
    ! KmMax が設定されていない場合. 
    ! 安定性解析では, dt / l**2 < 0.5 を満たす必要がある. ここでは 0.1 にしておいた. 
    !
    if (KmMax == 0.0d0) then 
       KmMax = 0.1 * (MixLen ** 2.0d0) / (DelTimeLong * 2.0d0) !LeapFrog
    end if

    !-------------------------------------------------------------------
    ! tendency の出力
    !
    call turbulence_kw1978_output

    !-------------------------------------------------------------------
    ! 配列の用意
    !
    allocate( xyr_DPTempBZDz(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyrf_DQMixBZDz(imin:imax, jmin:jmax, kmin:kmax, ncmax) )

    xyr_DPTempBZDz = xyr_dz_xyz( xyz_PTempBZ )
    do f = 1, ncmax
      xyrf_DQMixBZDz(:,:,:,f) = xyr_dz_xyz( xyzf_QMixBZ(:,:,:,f) )
    end do

    !-------------------------------------------------------------------
    ! Output
    !
    call MessageNotify( "M", &
      & module_name, "Cm = %f", d=(/Cm/))
    call MessageNotify( "M", &
      & module_name, "KmMax = %f", d=(/KmMax/))
    call MessageNotify( "M", &
      & module_name, "MixLen = %f", d=(/MixLen/))
    call MessageNotify( "M", &
      & module_name, "FlagArare4 = %b", l=(/ FlagArare4 /))
    call MessageNotify( "M", &
      & module_name, "FactorDExnerDtTurb= %f", d=(/ FactorDExnerDtTurb /))

  end subroutine turbulence_kw1978_init

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine turbulence_KW1978_forcing(       &
    & pyz_VelXBl,  xqz_VelYBl,  xyr_VelZBl,   &
    & xyz_PTempBl, xyz_ExnerBl, xyzf_QMixBl,  &
    & xyz_KmBl,    xyz_KhBl,    xyz_CDensBl,  &
    & pyz_DVelXDt, xqz_DVelYDt, xyr_DVelZDt,  &
    & xyz_DPTempDt,xyz_DExnerDt,xyzf_DQMixDt, &
    & xyz_DKmDt,   xyz_DCDensDt               &
    & )

    !モジュール読み込み 
    !
    use dc_types, only: DP, STRING    
    use gtool_historyauto, only: HistoryAutoPut
    use composition, only: SpcWetSymbol
    use timeset, only:  TimeN
    use gridset, only:  imin,           &! x 方向の配列の下限
      &                 imax,           &! x 方向の配列の上限
      &                 jmin,           &! y 方向の配列の下限
      &                 jmax,           &! y 方向の配列の上限
      &                 kmin,           &! z 方向の配列の下限
      &                 kmax,           &! z 方向の配列の上限
      &                 nx,ny,nz,ncmax
    use basicset, only: xyz_PTempBZ,     &!基本場の温位
      &                 xyzf_QMixBZ,        &!基本場の混合比 
      &                 xyz_ExnerBZ,        &!基本場のエクスナー関数
      &                 xyz_DensBZ,         &!基本場の密度
      &                 pyz_DensBZ,         &!基本場の密度
      &                 xqz_DensBZ,         &!基本場の密度
      &                 xyr_DensBZ           !基本場の密度
    use constants,only: Grav,          &
      &                 MolWtDry,      & 
      &                 CpDry
    use average,  only: xyz_pyz, xyr_pyr, xqz_pqz, &
      &                 pyz_xyz, pyr_xyr, pqz_xqz, &
      &                 xyz_xqz, pyz_pqz, xyr_xqr, &
      &                 xqz_xyz, pqz_pyz, xqr_xyr, &
      &                 xyz_xyr, pyz_pyr, xqz_xqr, &
      &                 xyr_xyz, pyr_pyz, xqr_xqz, &
      &                 pqz_xyz, pyr_xyz, xqr_xyz, &
      &                 xyz_pqz, xyz_pyr, xyz_xqr
    use differentiate_center2,  &
      &           only: xyz_dx_pyz, xyz_dy_xqz, xyz_dz_xyr, &
      &                 pyz_dx_xyz, xqz_dy_xyz, xyr_dz_xyz, &
      &                 pqz_dx_xqz, pqz_dy_pyz, pyz_dy_pqz, &
      &                 pyr_dx_xyr, pyz_dz_pyr, pyr_dz_pyz, &
      &                 pyr_dz_pyz, pyr_dx_xyr, xyr_dx_pyr, &
      &                 xqr_dz_xqz, xqr_dy_xyr, xyr_dy_xqr, &
      &                 pqz_dy_pyz, pqz_dx_xqz, xqz_dx_pqz, &
      &                 xqr_dy_xyr, xqr_dz_xqz, xqz_dz_xqr  
    use DExnerDt, only: xyz_DExnerDt_xyz

    ! 暗黙の型宣言禁止
    ! 
    implicit none

    ! 入出力変数
    !
    real(DP),intent(in) :: pyz_VelXBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !水平速度
    real(DP),intent(in) :: xqz_VelYBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !水平速度
    real(DP),intent(in) :: xyr_VelZBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !鉛直速度
    real(DP),intent(in) :: xyz_PTempBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !温位
    real(DP),intent(in) :: xyz_ExnerBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !無次元圧力
    real(DP),intent(in) :: xyzf_QMixBl(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                                    !凝縮成分の混合比
    real(DP),intent(in) :: xyz_KmBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !乱流拡散係数
    real(DP),intent(in) :: xyz_KhBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !乱流拡散係数
    real(DP),intent(in) :: xyz_CDensBl(imin:imax,jmin:jmax,kmin:kmax)

    real(DP),intent(inout):: pyz_DVelXDt(imin:imax,jmin:jmax,kmin:kmax)
                                                    !スカラー量の水平乱流拡散
    real(DP),intent(inout):: xqz_DVelYDt(imin:imax,jmin:jmax,kmin:kmax)
                                                    !スカラー量の水平乱流拡散
    real(DP),intent(inout):: xyr_DVelZDt(imin:imax,jmin:jmax,kmin:kmax)
                                                    !スカラー量の水平乱流拡散
    real(DP),intent(inout):: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)

    real(DP),intent(inout):: xyz_DExnerDt(imin:imax,jmin:jmax,kmin:kmax)

    real(DP),intent(inout):: xyzf_DQMixDt(imin:imax,jmin:jmax,kmin:kmax, ncmax)

    real(DP),intent(inout):: xyz_DKmDt(imin:imax,jmin:jmax,kmin:kmax)

    real(DP),intent(inout) :: xyz_DCDensDt(imin:imax,jmin:jmax,kmin:kmax)

    ! 作業変数
    !
    real(DP)            :: xyz_Buoy(imin:imax,jmin:jmax,kmin:kmax)
                                                    !渦粘性係数の
    real(DP)            :: xyz_BuoyT(imin:imax,jmin:jmax,kmin:kmax)
                                                    !渦粘性係数の
    real(DP)            :: xyz_BuoyM(imin:imax,jmin:jmax,kmin:kmax)
                                                    !渦粘性係数の
    real(DP)            :: xyz_Shear(imin:imax,jmin:jmax,kmin:kmax)
                                                    !渦粘性係数の
    real(DP)            :: xyz_Diff(imin:imax,jmin:jmax,kmin:kmax)
                                                    !渦粘性係数の
    real(DP)            :: xyz_Disp(imin:imax,jmin:jmax,kmin:kmax)
                                                    !乱流エネルギーの消散
    real(DP)            :: xyz_DispPI(imin:imax,jmin:jmax,kmin:kmax)
                                                    !乱流エネルギーの消散
    real(DP)            :: xyz_DispHeat(imin:imax,jmin:jmax,kmin:kmax)
                                                    !乱流エネルギーの消散
    real(DP)            :: xyz_Turb(imin:imax,jmin:jmax,kmin:kmax)
                                                    !
    real(DP)            :: xyzf_Turb(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                                    !スカラー量の水平乱流拡散
    real(DP)            :: pyz_Turb(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: xqz_Turb(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: xyr_Turb(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: pyz_DVelXDt0(imin:imax,jmin:jmax,kmin:kmax)
                                                    !スカラー量の水平乱流拡散
    real(DP)            :: xqz_DVelYDt0(imin:imax,jmin:jmax,kmin:kmax)
                                                    !スカラー量の水平乱流拡散
    real(DP)            :: xyr_DVelZDt0(imin:imax,jmin:jmax,kmin:kmax)
                                                    !スカラー量の水平乱流拡散
    real(DP)            :: xyz_DPTempDt0(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: xyz_DExnerDt0(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: xyzf_DQMixDt0(imin:imax,jmin:jmax,kmin:kmax, ncmax)

    real(DP)            :: xyz_DKmDt0(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: xyz_DCDensDt0(imin:imax,jmin:jmax,kmin:kmax)

    integer             :: f


    !----------------------------------
    ! 初期化
    !
    pyz_DVelXDt0 = pyz_DVelXDt
    xqz_DVelYDt0 = xqz_DVelYDt
    xyr_DVelZDt0 = xyr_DVelZDt
    xyz_DPTempDt0 = xyz_DPTempDt
    xyz_DExnerDt0 = xyz_DExnerDt
    xyzf_DQMixDt0 = xyzf_DQMixDt
    xyz_DKmDt0    = xyz_DKmDt
    xyz_DCDensDt0 = xyz_DCDensDt
    
    !----------------------------------
    ! 拡散係数の時間発展 (エネルギー方程式を Km の式に変形したもの)
    !

    ! Buoyancy term
    !
    xyz_Buoy = xyz_BuoyMoistKm(xyz_PTempBl, xyz_ExnerBl, xyzf_QMixBl) 

    xyz_BuoyT =                                                  &
      &  - 3.0d0 * Grav * ( Cm * Cm  * MixLen * MixLen )         &
      &    * xyz_xyr(                                            &
      &        xyr_dz_xyz( xyz_PTempBl ) + xyr_DPTempBZDz        &
      &      ) / ( 2.0d0 * xyz_PTempBZ )
    
    xyz_BuoyM = xyz_Buoy - xyz_BuoyT

    xyz_Shear = &
      &   ( ( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ) )                &
      & * (                                                        &
      &      ( xyz_dx_pyz( pyz_VelXBl ) ) ** 2.0d0                 &
      &    + ( xyz_dy_xqz( xqz_VelYBl ) ) ** 2.0d0                 &
      &    + ( xyz_dz_xyr( xyr_VelZBl ) ) ** 2.0d0                 &
      &    + 5.0d-1                                                &
      &      * (                                                   &
      &          (                                                 &
      &              xyz_pyr( pyr_dz_pyz( pyz_VelXBl ) )           &
      &            + xyz_pyr( pyr_dx_xyr( xyr_VelZBl ) )           &
      &           ) ** 2.0d0                                       &
      &        + (                                                 &
      &              xyz_xqr( xqr_dy_xyr( xyr_VelZBl ) )           &
      &            + xyz_xqr( xqr_dz_xqz( xqz_VelYBl ) )           &
      &           ) ** 2.0d0                                       &
      &        + (                                                 &
      &              xyz_pqz( pqz_dx_xqz( xqz_VelYBl ) )           &
      &            + xyz_pqz( pqz_dy_pyz( pyz_VelXBl ) )           &
      &           ) ** 2.0d0                                       &
      &        )                                                   &
      &   )                                                        &
      & - xyz_KmBl * (  xyz_dx_pyz( pyz_VelXBl )                   &
      &               + xyz_dy_xqz( xqz_VelYBl )                   &
      &               + xyz_dz_xyr( xyr_VelZBl ) ) / 3.0d0

    xyz_Diff =                                           &
      & + (                                              &
      &    + xyz_dx_pyz(pyz_dx_xyz(xyz_KmBl ** 2.0d0))   &
      &    + xyz_dy_xqz(xqz_dy_xyz(xyz_KmBl ** 2.0d0))   &
      &    + xyz_dz_xyr(xyr_dz_xyz(xyz_KmBl ** 2.0d0))   &
      &   ) * 5.0d-1                                     &
      & + (                                              &
      &    + (xyz_pyz(pyz_dx_xyz(xyz_KmBl))) ** 2.0d0    &
      &    + (xyz_xqz(xqz_dy_xyz(xyz_KmBl))) ** 2.0d0    &
      &    + (xyz_xyr(xyr_dz_xyz(xyz_KmBl))) ** 2.0d0    &
      &   )

    ! t - \Delta t で評価
    !
    xyz_Disp = - (xyz_KmBl ** 2.0d0) * 5.0d-1 / (MixLen ** 2.0d0)

    ! tendency
    !
    xyz_DKmDt = ( xyz_Buoy + xyz_Shear + xyz_Diff + xyz_Disp ) + xyz_DKmDt0 

    call HistoryAutoPut(TimeN, 'DKmDtBuoy',  xyz_Buoy(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DKmDtBuoyT', xyz_BuoyT(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DKmDtBuoyM', xyz_BuoyM(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DKmDtShear', xyz_Shear(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DKmDtTurb',  xyz_Diff(1:nx,1:ny,1:nz)) 
    call HistoryAutoPut(TimeN, 'DKmDtDisp',  xyz_Disp(1:nx,1:ny,1:nz))    
    
    !--------------------------------
    ! 雲密度の tendency
    !
    if (FlagArare4) then 
    xyz_Turb =                                                           &
      &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyz_CDensBl ) )  &
      & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyz_CDensBl ) )  &
      & + xyz_dz_xyr( xyr_xyz( xyz_KhBl ) * xyr_dz_xyz( xyz_CDensBl ) )
    else
    xyz_Turb =                                                           &
      &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyz_CDensBl ) )  &
      & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyz_CDensBl ) )  &
      & + xyz_dz_xyr(                                                    &
      &     xyr_xyz( xyz_DensBZ * xyz_KhBl ) * xyr_dz_xyz( xyz_CDensBl ) &
      &   ) / xyz_DensBZ
    end if

    xyz_DCDensDt = xyz_DCDensDt0 + xyz_Turb

    call HistoryAutoPut(TimeN, 'DCDensDtTurb', xyz_Turb(1:nx,1:ny,1:nz))

    !--------------------------------
    ! 温位の tendency
    !
    if (FlagArare4) then 
    xyz_Turb =                                                           &
      &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyz_PTempBl ) )  &
      & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyz_PTempBl ) )  &
      & + xyz_dz_xyr( xyr_xyz( xyz_KhBl ) * xyr_dz_xyz( xyz_PTempBl ) )  &
      & + xyz_dz_xyr( xyr_xyz( xyz_KhBl ) * xyr_DPTempBZDz )
    else
    xyz_Turb =                                                           &
      &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyz_PTempBl ) )  &
      & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyz_PTempBl ) )  &
      & + xyz_dz_xyr(                                                    &
      &     xyr_xyz( xyz_DensBZ * xyz_KhBl )                             &
      &     * ( xyr_dz_xyz( xyz_PTempBl ) + xyr_DPTempBZDz )             &
      &   )                                                              &
      &   / xyz_DensBZ
    end if
    
    xyz_DispHeat = (xyz_KmBl ** 3.0d0)                                   &
      & / (xyz_ExnerBZ * CpDry * (Cm ** 2.0d0) * (MixLen ** 4.0d0))

    xyz_DPTempDt = ( xyz_Turb + xyz_DispHeat ) + xyz_DPTempDt0 

    call HistoryAutoPut(TimeN, 'DPTempDtDisp', xyz_DispHeat(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'DPTempDtTurb', xyz_Turb(1:nx, 1:ny, 1:nz))

    !--------------------------------
    ! 混合比の tendency
    !
    do f = 1, ncmax    
      if (FlagArare4) then 
      xyzf_Turb(:,:,:,f) =                                                         &
        & + xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyzf_QMixBl(:,:,:,f) ) ) &
        & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyzf_QMixBl(:,:,:,f) ) ) &
        & + xyz_dz_xyr( xyr_xyz( xyz_KhBl ) * xyr_dz_xyz( xyzf_QMixBl(:,:,:,f) ) ) &
        & + xyz_dz_xyr( xyr_xyz( xyz_KhBl ) * xyrf_DQMixBZDz(:,:,:,f) )
      else
      xyzf_Turb(:,:,:,f) =                                                         &
        &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyzf_QMixBl(:,:,:,f) ) ) &
        & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyzf_QMixBl(:,:,:,f) ) ) &
        & + xyz_dz_xyr(                                                            &
        &     xyr_xyz( xyz_DensBZ * xyz_KhBl )                                     &
        &     * ( xyr_dz_xyz( xyzf_QMixBl(:,:,:,f) ) + xyrf_DQMixBZDz(:,:,:,f) )   &
        &   ) / xyz_DensBZ 
      end if 
      
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtTurb', xyzf_Turb(1:nx,1:ny,1:nz,f))
    end do

    xyzf_DQMixDt = xyzf_DQMixDt0 + xyzf_Turb
    
    !--------------------------------
    ! VelX の tendency
    !
    if (FlagArare4) then 
    pyz_Turb =                                                     &
      &   2.0d0 * pyz_dx_xyz( xyz_KmBl * xyz_dx_pyz( pyz_VelXBl ) )&
      & + pyz_dy_pqz(                                              &
      &       pqz_xyz( xyz_KmBl ) * pqz_dx_xqz( xqz_VelYBl )       &
      &     + pqz_xyz( xyz_KmBl ) * pqz_dy_pyz( pyz_VelXBl )       &
      &   )                                                        &
      & + pyz_dz_pyr(                                              &
      &       pyr_xyz( xyz_KmBl ) * pyr_dx_xyr( xyr_VelZBl )       &
      &     + pyr_xyz( xyz_KmBl ) * pyr_dz_pyz( pyz_VelXBl )       &
      &   )                                                        &
      & - 2.0d0 * pyz_dx_xyz( ( xyz_KmBl ** 2.0d0 ) )              &
      &   / ( 3.0d0 * ( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ) )
    else
    pyz_Turb =                                                             &
      &   2.0d0 * pyz_dx_xyz( xyz_KmBl * xyz_dx_pyz( pyz_VelXBl ) )        &
      & + pyz_dy_pqz(                                                      &
      &       pqz_xyz( xyz_KmBl ) * pqz_dx_xqz( xqz_VelYBl )               &
      &     + pqz_xyz( xyz_KmBl ) * pqz_dy_pyz( pyz_VelXBl )               &
      &   )                                                                &
      & + pyz_dz_pyr(                                                      &
      &       pyr_xyz( xyz_DensBZ * xyz_KmBl ) * pyr_dx_xyr( xyr_VelZBl )  &
      &     + pyr_xyz( xyz_DensBZ * xyz_KmBl ) * pyr_dz_pyz( pyz_VelXBl )  &
      &   ) / pyz_DensBZ                                                   &
      & - 2.0d0 * pyz_dx_xyz( ( xyz_KmBl ** 2.0d0 ) )                      &
      &   / ( 3.0d0 * ( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ) )
    end if

    pyz_DVelXDt = pyz_DVelXDt0 + pyz_Turb

    call HistoryAutoPut(TimeN, 'DVelXDtTurb', pyz_Turb(1:nx, 1:ny, 1:nz))


    !--------------------------------
    ! VelY の tendency
    !
    if (FlagArare4) then 
    xqz_Turb =                                                      &
      &   2.0d0 * xqz_dy_xyz( xyz_KmBl * xyz_dy_xqz( xqz_VelYBl ) ) &
      & + xqz_dx_pqz(                                               &
      &       pqz_xyz( xyz_KmBl ) * pqz_dy_pyz( pyz_VelXBl )        &
      &     + pqz_xyz( xyz_KmBl ) * pqz_dx_xqz( xqz_VelYBl )        &
      &   )                                                         &
      & + xqz_dz_xqr(                                               &
      &       xqr_xyz( xyz_KmBl ) * xqr_dy_xyr( xyr_VelZBl )        &
      &     + xqr_xyz( xyz_KmBl ) * xqr_dz_xqz( xqz_VelYBl )        &
      &   )                                                         &
      & - 2.0d0 * xqz_dy_xyz( ( xyz_KmBl ** 2.0d0 ) )               &
      &   / ( 3.0d0 * ( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ) )
    else
    xqz_Turb =                                                            &
      &   2.0d0 * xqz_dy_xyz( xyz_KmBl * xyz_dy_xqz( xqz_VelYBl ) )       &
      & + xqz_dx_pqz(                                                     &
      &       pqz_xyz( xyz_KmBl ) * pqz_dy_pyz( pyz_VelXBl )              &
      &     + pqz_xyz( xyz_KmBl ) * pqz_dx_xqz( xqz_VelYBl )              &
      &   )                                                               &
      & + xqz_dz_xqr(                                                     &
      &       xqr_xyz( xyz_DensBZ * xyz_KmBl ) * xqr_dy_xyr( xyr_VelZBl ) &
      &     + xqr_xyz( xyz_DensBZ * xyz_KmBl ) * xqr_dz_xqz( xqz_VelYBl ) &
      &   ) / xqz_DensBZ                                                  &
      & - 2.0d0 * xqz_dy_xyz( ( xyz_KmBl ** 2.0d0 ) )                     &
      &   / ( 3.0d0 * ( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ) )
    end if 
    
    xqz_DVelYDt = xqz_DVelYDt0 + xqz_Turb

    call HistoryAutoPut(TimeN, 'DVelYDtTurb', xqz_Turb(1:nx, 1:ny, 1:nz))


    !--------------------------------
    ! VelZ の tendency
    !
    if (FlagArare4) then 
    xyr_Turb =                                                       &
      & + 2.0d0 * xyr_dz_xyz( xyz_KmBl * xyz_dz_xyr( xyr_VelZBl ) )  &
      & + xyr_dx_pyr(                                                &
      &      pyr_xyz( xyz_KmBl ) * pyr_dz_pyz( pyz_VelXBl )          &
      &    + pyr_xyz( xyz_KmBl ) * pyr_dx_xyr( xyr_VelZBl )          &
      &   )                                                          &
      & + xyr_dy_xqr(                                                &
      &      xqr_xyz( xyz_KmBl ) * xqr_dz_xqz( xqz_VelYBl )          &
      &    + xqr_xyz( xyz_KmBl ) * xqr_dy_xyr( xyr_VelZBl )          &
      &   )                                                          & 
      & - 2.0d0 * xyr_dz_xyz(  xyz_KmBl ** 2.0d0 )                   &
      &   / ( 3.0d0 * ( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ) )
    else
    xyr_Turb =                                                       &
      & + 2.0d0                                                      &
      &   * xyr_dz_xyz(                                              &
      &       xyz_DensBZ * xyz_KmBl * xyz_dz_xyr( xyr_VelZBl )       &
      &     ) / xyr_DensBZ                                           &
      & + xyr_dx_pyr(                                                &
      &      pyr_xyz( xyz_KmBl ) * pyr_dz_pyz( pyz_VelXBl )          &
      &    + pyr_xyz( xyz_KmBl ) * pyr_dx_xyr( xyr_VelZBl )          &
      &   )                                                          &
      & + xyr_dy_xqr(                                                &
      &      xqr_xyz( xyz_KmBl ) * xqr_dz_xqz( xqz_VelYBl )          &
      &    + xqr_xyz( xyz_KmBl ) * xqr_dy_xyr( xyr_VelZBl )          &
      &   )                                                          & 
      & - 2.0d0 * xyr_dz_xyz( xyz_DensBZ * xyz_KmBl * xyz_KmBl )     &
      &   / ( 3.0d0 * ( Cm * Cm * MixLen * MixLen ) )                &
      &   / xyr_DensBZ 
    end if

    xyr_DVelZDt = xyr_DVelZDt0 + xyr_Turb

    call HistoryAutoPut(TimeN, 'DVelZDtTurb', xyr_Turb(1:nx, 1:ny, 1:nz))

    !--------------------
    ! Exner function
    !
    xyz_DispPI = xyz_DExnerDt_xyz( xyz_DispHeat ) * FactorDExnerDtTurb
    xyz_DExnerDt = xyz_DExnerDt0 + xyz_DispPI

    call HistoryAutoPut(TimeN, 'DExnerDtDisp', xyz_DispPI(1:nx, 1:ny, 1:nz))


  contains

  function xyz_BuoyMoistKm(xyz_PTemp, xyz_Exner, xyzf_QMix)
    !
    ! 浮力項の計算
    !

    !モジュール呼び出し
    !
    use dc_types, only: DP, STRING    
    use composition, only: MolWtWet,        &
      &                 SpcWetID,           &
      &                 CondNum,            &!凝結過程の数
      &                 IdxCG,              &!凝結過程(蒸気)の配列添え字
      &                 IdxCC,              &!凝結過程(雲)の配列添え字
      &                 GasNum,             &!気体の数
      &                 IdxG                 !気体の配列添え字
    use ChemCalc, only: xyz_LatentHeat       !潜熱
    !    &              ReactHeatNH4SH       !NH4SH の反応熱
    use gtool_historyauto,                  &
      &           only: HistoryAutoPut
    use basicset, only: xyz_PTempBZ,        &!基本場の温位
      &                 xyz_QMixBZ,         &!基本場の混合比
      &                 xyz_QMixBZPerMolWt, &!基本場の混合比
      &                 xyz_EffMolWtBZ,     &!基本場の混合比
      &                 xyz_ExnerBZ          !基本場のエクスナー関数
    use average,  only: xyz_xyr
    use differentiate_center2,              &
      &           only: xyr_dz_xyz

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in) :: xyz_PTemp(imin:imax,jmin:jmax,kmin:kmax)
                                               !温位
    real(DP), intent(in) :: xyz_Exner(imin:imax,jmin:jmax,kmin:kmax)
                                               !無次元圧力
    real(DP), intent(in) :: xyzf_QMix(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                               !凝縮成分の混合比
    real(DP) :: xyzf_QMixAll(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                               !凝縮成分の混合比
    real(DP) :: xyzf_QMixAll2(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                               !凝縮成分の混合比
    real(DP) :: xyz_BuoyMoistKm(imin:imax,jmin:jmax,kmin:kmax)
                                               !
    real(DP) :: xyzf_LatentHeat(imin:imax,jmin:jmax,kmin:kmax, GasNum)
                                               !潜熱
    real(DP) :: xyz_TempAll(imin:imax,jmin:jmax,kmin:kmax)
                                               !温度
    real(DP) :: xyz_EffHeat(imin:imax,jmin:jmax,kmin:kmax)
                                               !
    real(DP) :: xyz_EffPTemp(imin:imax,jmin:jmax,kmin:kmax)    
                                               !浮力に対する温度差の寄与
    real(DP) :: xyz_EffMolWt(imin:imax,jmin:jmax,kmin:kmax)    
                                               !浮力に対する分子量差の寄与
    real(DP) :: xyzf_QMixPerMolWt(imin:imax,jmin:jmax,kmin:kmax, GasNum)
                                               !混合比/分子量
    integer              :: s
    
    
    !温度, 圧力, 混合比の全量を求める
    !擾乱成分と平均成分の足し算
    xyz_TempAll     = ( xyz_PTemp + xyz_PTempBZ ) * ( xyz_Exner + xyz_ExnerBZ )
    xyzf_QMixAll    = xyzf_QMixBZ + xyzf_QMix
    xyzf_LatentHeat = 0.0d0
    
    !作業配列の初期化. 気体のみ利用
    do s = 1, GasNum
      xyzf_QMixPerMolWt(:,:,:,s) = xyzf_QMix(:,:,:,IdxG(s)) / MolWtWet(IdxG(s))
    end do
    
    !温度の効果
    xyz_EffPTemp = xyz_PTemp / xyz_PTempBZ 
    
    !分子量効果 + 引きづりの効果
    xyz_EffMolWt =                                      &
      & + sum(xyzf_QMixPerMolWt, 4)                     &
      &    / ( 1.0d0 / MolWtDry + xyz_QMixBZPerMolWt )  &
      & - sum(xyzf_QMix, 4) / ( 1.0d0 + xyz_QMixBZ )
    
    !蒸気が蒸発する場合の潜熱を計算
    !  分子量の部分はいつでも効くが潜熱は飽和していないと効かないので, 
    !  雲の混合比がゼロの時には, 潜熱の寄与はゼロとなるように調節している
    !
    xyzf_QMixAll2 = xyzf_QMixAll
!    xzya_QMixAll2(:,:,:,IdxNH3) = &
!      &  xyzf_QMixAll(:,:,:,IdxNH3) - xyzf_QMixAll(:,:,:,IdxH2S) 

    do s = 1, CondNum
      xyzf_LatentHeat(:,:,:,s) =                                                   &
        & xyz_LatentHeat( SpcWetID(IdxCC(s)), xyz_TempAll )                        &
        &  * xyzf_QMixAll2(:,:,:,IdxCG(s))                                         &
        &  * ( 5.0d-1 + sign( 5.0d-1, (xyzf_QMixAll2(:,:,:,IdxCC(s)) - 1.0d-4) ) )
    end do

    xyz_EffHeat = ( sum( xyzf_LatentHeat, 4 ) * xyz_EffMolWtBZ &
!      &            + ReactHeatNH4SH * &
!      &            (xyzf_QMixAll(:,:,:,IdxH2S) + xyzf_QMixAll(:,:,:,IdxH2S)) &
      &          ) / ( CpDry * xyz_ExnerBZ ) 

    
    !乱流拡散係数の時間発展式の浮力項を決める
    xyz_BuoyMoistKm = &
      &  - 3.0d0 * Grav * (( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ))   &
      &    * xyz_xyr(                                               &
      &        xyr_dz_xyz(                                          &
      &           xyz_EffHeat                                       &
      &         + xyz_PTempBZ / xyz_EffMolWtBZ                      &
      &           * ( 1.0d0 + xyz_EffPTemp + xyz_EffMolWt )         &
      &         )                                                   &
      &      )                                                      &
      &    / ( 2.0d0 * xyz_PTempBZ / xyz_EffMolWtBZ)                   

  end function xyz_BuoyMoistKm

  end subroutine turbulence_KW1978_forcing

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine turbulence_kw1978_output
    !
    ! tendency の出力のための設定
    !

    ! モジュール呼び出し
    !
    use gtool_historyauto, only : HistoryAutoAddVariable
    use gridset,           only : ncmax
    use composition,       only : SpcWetSymbol

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 作業変数
    !
    integer :: l

    !-----------------------------------------------------
    ! tendency 出力のための設定
    !
    call HistoryAutoAddVariable(  &
      & varname='DVelXDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of velocity (x)', &
      & units='m.s-2',    &
      & xtype='double')
    
    call HistoryAutoAddVariable(  &
      & varname='DVelYDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of velocity (y)', &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of velocity (z)', &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DPTempDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of potential temperature', &
      & units='K.s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DPTempDtDisp',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Dissipation term of potential temperature', &
      & units='K.s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DExnerDtDisp',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Dissipation term of exner function', &
      & units='s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DCDensDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of cloud density', &
      & units='kg.m-3.s-1',    &
      & xtype='double')

    do l = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtTurb', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='Turbulence term of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio',  &
        & units='kg.kg-1.s-1',    &
        & xtype='double')
    end do

    call HistoryAutoAddVariable(  &
      & varname='DKmDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Diffusion term of Km', &
      & units='s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DKmDtBuoy',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Buoyancy term of Km', &
      & units='s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DKmDtBuoyT',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Buoyancy term of Km (temperature)', &
      & units='s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DKmDtBuoyM',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Buoyancy term of Km (composition)', &
      & units='s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DKmDtShear',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Shear term of Km', &
      & units='s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DKmDtDisp',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Dissipation term of Km', &
      & units='s-1',    &
      & xtype='double')

  end subroutine turbulence_kw1978_output

end module Turbulence_kw1978_v1
