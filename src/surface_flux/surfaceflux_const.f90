!= Module HeatFluxConst
!
! Authors::   ODAKA Masatsugu
! Version::   $Id: surfaceflux_const.f90,v 1.5 2014/03/04 04:49:42 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2012. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Surfaceflux_const
  !
  ! 上部, 下部境界からの一定の熱・運動量・物質フラックスを与えた場合による
  ! 運動量, 温度, 凝結成分の変化率を計算するモジュール. 
  !

  !モジュール読み込み
  use dc_types, only: DP, STRING

  !暗黙の型宣言禁止
  implicit none

  !属性の指定
  private

  !関数を public に設定
  public surfaceflux_const_init
  public surfaceflux_const_forcing

  !変数定義
  real(DP), save :: SfcXMomFluxBtm = 0.0d0   ! X 方向の運動量フラックス (下部境界) 
  real(DP), save :: SfcXMomFluxTop = 0.0d0   ! X 方向の運動量フラックス (上部境界) 
  real(DP), save :: SfcYMomFluxBtm = 0.0d0   ! Y 方向の運動量フラックス (下部境界) 
  real(DP), save :: SfcYMomFluxTop = 0.0d0   ! Y 方向の運動量フラックス (上部境界) 
  real(DP), save :: SfcHeatFluxBtm = 0.0d0   ! 熱フラックス (下部境界) 
  real(DP), save :: SfcHeatFluxTop = 0.0d0   ! 熱フラックス (上部境界) 
  real(DP), save :: SfcQmixFluxBtm = 0.0d0   ! 物質フラックス (下部境界) 
  real(DP), save :: SfcQmixFluxTop = 0.0d0   ! 物質フラックス (上部境界) 

  character(STRING), parameter:: module_name = 'surfaceflux_const'
                              ! モジュールの名称.
                              ! Module name
  real(DP), save :: FactorDExnerDtSurf = 1.0d0
 

contains
!!!------------------------------------------------------------------------!!!
  subroutine Surfaceflux_const_init
    !
    !NAMELIST から必要な情報を読み取り, 出力変数の設定を行う. 
    !

    !モジュール読み込み
    use dc_iounit,         only : FileOpen
    use dc_message,        only : MessageNotify
    use gtool_historyauto, only : HistoryAutoAddVariable    
    use gridset,           only : ncmax
    use composition,       only : SpcWetSymbol
    use namelist_util,     only : namelist_filename

    !暗黙の型宣言禁止
    implicit none

    !内部変数
    integer    :: l, unit

    !---------------------------------------------------------------
    ! NAMELIST から情報を取得
    !
    NAMELIST /surfaceflux_const_nml/ &
      & SfcXMomFluxBtm, SfcXMomFluxTop, SfcYMomFluxBtm, SfcYMomFluxTop, &
      & SfcHeatFluxBtm, SfcHeatFluxTop, SfcQmixFluxBtm, SfcQmixFluxTop, &
      & FactorDExnerDtSurf 

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=surfaceflux_const_nml)
    close(unit)  

    call MessageNotify( "M", module_name, "SfcXMomFluxBtm = %f", &
      &                  d=(/SfcXMomFluxBtm/) )
    call MessageNotify( "M", module_name, "SfcXMomFluxTop = %f", &
      &                  d=(/SfcXMomFluxTop/) )
    call MessageNotify( "M", module_name, "SfcYMomFluxBtm = %f", &
      &                  d=(/SfcYMomFluxBtm/) )
    call MessageNotify( "M", module_name, "SfcYMomFluxTop = %f", &
      &                  d=(/SfcYMomFluxTop/) )
    call MessageNotify( "M", module_name, "SfcHeatFluxBtm = %f", &
      &                  d=(/SfcHeatFluxBtm/) )
    call MessageNotify( "M", module_name, "SfcHeatFluxTop = %f", &
      &                  d=(/SfcHeatFluxTop/) )
    call MessageNotify( "M", module_name, "SfcQmixFluxBtm = %f", &
      &                  d=(/SfcQmixFluxBtm/) )
    call MessageNotify( "M", module_name, "SfcQmixFluxTop = %f", &
      &                  d=(/SfcQmixFluxTop/) )
    call MessageNotify( 'M', module_name, "FactorDExnerDtSurf = %f", &
      &                  d=(/FactorDExnerDtSurf/) )

    !---------------------------------------------------------------
    ! ファイル出力の定義
    !
    call HistoryAutoAddVariable(      &
      & varname='PTempSfcFlux',       &
      & dims=(/'x','y','z','t'/),         &
      & longname='surface potential temperature flux (heat flux divided by density and specific heat)', &
      & units='K.m.s-1',             &
      & xtype='float')

    call HistoryAutoAddVariable(      &
      & varname='ExnerSfcFlux',       &
      & dims=(/'x','y','z','t'/),         &
      & longname='surface exner function flux (heat flux divided by density and specific heat)', &
      & units='s-1',             &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='VelXSfcFlux',    &
      & dims=(/'x','y','z','t'/),     &
      & longname='surface flux of x-component of velocity (momentum flux divided by density)', &
      & units='m2.s-2',           &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='VelYSfcFlux',    &
      & dims=(/'x','y','z','t'/),     &
      & longname='surface flux of y-component of velocity (momentum flux divided by density)', &
      & units='m2.s-2',           &
      & xtype='float')

    do l = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname=trim(SpcWetSymbol(l))//'SfcFlux', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='surface flux of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio (mass flux divided by density)',  &
        & units='m.s-1',    &
        & xtype='float')
    end do


    call HistoryAutoAddVariable(  &
      & varname='DPTempDtSfc',         &
      & dims=(/'x','y','z','t'/), &
      & longname='potential temperature tendency by surface flux', &
      & units='K.s-1',            &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='DExnerDtSfc',         &
      & dims=(/'x','y','z','t'/), &
      & longname='exner function tendency by surface flux', &
      & units='K.s-1',            &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='DVelXDtSfc',         &
      & dims=(/'x','y','z','t'/), &
      & longname='x-component velocity tendency by surface flux', &
      & units='m.s-2',            &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='DVelYDtSfc',         &
      & dims=(/'x','y','z','t'/), &
      & longname='y-component velocity tendency by surface flux', &
      & units='m.s-2',            &
      & xtype='float')

    do l = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtSfc', & 
        & dims=(/'x','y','z','t'/),     &
        & longname=trim(SpcWetSymbol(l))//' mixing ratio tendency by surface flux',  &
        & units='s-1',    &
        & xtype='float')
    end do


    call HistoryAutoAddVariable(       &
      & varname='SfcHeatFlux',         &
      & dims=(/'x','y','z','t'/),          &
      & longname='surface heat flux',  &
      & units='W.m-2',                 &
      & xtype='float')

    call HistoryAutoAddVariable(                      &
      & varname='SfcXMomFlux',                        &
      & dims=(/'x','y','z','t'/),                         &
      & longname='surface x-component momentum flux', &
      & units='kg.m-2.s-1',                           &
      & xtype='float')

    call HistoryAutoAddVariable(                      &
      & varname='SfcYMomFlux',                        &
      & dims=(/'x','y','z','t'/),                         &
      & longname='surface y-component momentum flux', &
      & units='kg.m-2.s-1',                           &
      & xtype='float')

    do l = 1, ncmax
      call HistoryAutoAddVariable(                               &
        & varname=trim(SpcWetSymbol(l))//'_SfcMassFlux',         &
        & dims=(/'x','y','z','t'/),                              &
        & longname=trim(SpcWetSymbol(l))//' surface mass flux',  &
        & units='kg.m-2.s-1',                                    &
        & xtype='float')
    end do

  end subroutine Surfaceflux_Const_init


!!!------------------------------------------------------------------------!!!
  subroutine Surfaceflux_Const_forcing( &
    &   pyz_DVelXDt, xqz_DVelYDt, xyz_DPTempDt, xyz_DExnerDt, xyzf_DQMixDt &
    & )
    ! 

    !モジュール読み込み
    use dc_types,          only : DP
    use gtool_historyauto, only : HistoryAutoPut
    use gridset,           only : imin,         & !x 方向の配列の下限
      &                           imax,         & !x 方向の配列の上限
      &                           jmin,         & !y 方向の配列の下限
      &                           jmax,         & !y 方向の配列の上限
      &                           kmin,         & !z 方向の配列の下限
      &                           kmax,         & !z 方向の配列の上限
      &                           nx, ny, nz, ncmax
    use axesset,           only : xyz_dz          !z 方向の格子点間隔
!    use xyz_base_module,   only : xyz_pyz,  &
!      &                           xyz_xqz,  &
!      &                           pyz_xyz,  &
!      &                           xqz_xyz
    use basicset,          only : xyz_DensBZ     !基本場の密度
    use constants,         only : CpDry
    use composition,       only : SpcWetSymbol
    use timeset,           only : TimeN
    use DExnerDt,          only : xyz_DExnerDt_xyz_xyzf

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(inout):: pyz_DVelXDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout):: xqz_DVelYDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout):: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout):: xyz_DExnerDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout):: xyzf_DQMixDt(imin:imax,jmin:jmax,kmin:kmax, ncmax)

    real(DP) :: pyr_VelXflux (imin:imax,jmin:jmax,kmin:kmax) !速度フラックス
    real(DP) :: xqr_VelYflux (imin:imax,jmin:jmax,kmin:kmax) !速度フラックス
    real(DP) :: xyr_PTempFlux(imin:imax,jmin:jmax,kmin:kmax) !温度フラックス
    real(DP) :: xyr_ExnerFlux(imin:imax,jmin:jmax,kmin:kmax) !Exner 関数強制項
    real(DP) :: xyrf_QMixFlux(imin:imax,jmin:jmax,kmin:kmax,ncmax)
                                                             !物質フラックス
    real(DP) :: pyr_SfcXMomFlux(imin:imax,jmin:jmax,kmin:kmax) !運動量フラックス
    real(DP) :: xqr_SfcYMomFlux(imin:imax,jmin:jmax,kmin:kmax) !運動量フラックス
    real(DP) :: xyr_SfcHeatFlux(imin:imax,jmin:jmax,kmin:kmax) !熱フラックス

    real(DP) :: pyz_DVelXDtFlux (imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xqz_DVelYDtFlux (imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_DPTempDtFlux(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_DExnerDtFlux(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyzf_DQMixDtFlux(imin:imax,jmin:jmax,kmin:kmax, ncmax)

    integer  :: s             !ループ変数

!-----

    ! 初期化
    xyr_SfcHeatFlux = 0.0d0
    xyr_SfcHeatFlux = 0.0d0
    
    pyr_SfcXMomFlux = 0.0d0
    pyr_SfcXMomFlux = 0.0d0

    xqr_SfcYMomFlux = 0.0d0
    xqr_SfcYMomFlux = 0.0d0

    xyrf_QMixFlux = 0.0d0
    xyrf_QMixFlux = 0.0d0
    

    ! NAMLIST で設定した値を代入
    !
    xyr_SfcHeatFlux(:,:,1)  = SfcHeatFluxBtm
    xyr_SfcHeatFlux(:,:,nz) = SfcHeatFluxTop
    
    pyr_SfcXMomFlux(:,:,1)  = SfcXMomFluxBtm
    pyr_SfcXMomFlux(:,:,nz) = SfcXMomFluxTop

    xqr_SfcYMomFlux(:,:,1)  = SfcYMomFluxBtm
    xqr_SfcYMomFlux(:,:,nz) = SfcYMomFluxTop

    xyrf_QMixFlux(:,:,1,:)  = SfcQmixFluxBtm
    xyrf_QMixFlux(:,:,nz,:) = SfcQmixFluxTop

    ! 速度フラックス, 温度フラックスへの変換
    ! 
    pyr_VelXFlux  = pyr_SfcXMomFlux / xyz_DensBz
    xqr_VelYFlux  = xqr_SfcYMomFlux / xyz_DensBz
    xyr_PtempFlux = xyr_SfcHeatFlux / xyz_DensBz / CpDry

    ! Exner 関数予報式の強制項の計算
    !
    xyr_ExnerFlux = xyz_DExnerDt_xyz_xyzf(xyr_PTempFlux, xyrf_QMixFlux) * FactorDExnerDtSurf 

    ! 時間変化率の計算
    !
    xyz_DPTempDtFlux = xyr_PTempFlux / xyz_dz
    xyz_DExnerDtFlux = xyr_ExnerFlux / xyz_dz
    do s = 1, ncmax
      xyzf_DQMixDtFlux(:,:,:,s) = xyrf_QMixFlux(:,:,:,s) / xyz_dz
    end do
    pyz_DVelXDtFlux  = pyr_VelXFlux / xyz_dz
    xqz_DVelYDtFlux  = xqr_VelYFlux / xyz_dz

    ! 仮引数配列への加算
    !
    xyz_DPTempDt = xyz_DPTempDt + xyz_DPTempDtFlux
    xyz_DExnerDt = xyz_DExnerDt + xyz_DExnerDtFlux
    do s = 1, ncmax
      xyzf_DQMixDt(:,:,:,s) = xyzf_DQMixDt(:,:,:,s) + xyzf_DQMixDtFlux(:,:,:,s)
    end do
    pyz_DVelXDt = pyz_DVelXDt + pyz_DVelXDtFlux
    xqz_DVelYDt = xqz_DVelYDt + xqz_DVelYDtFlux

    ! 出力
    ! 
    call HistoryAutoPut(TimeN, 'DPTempDtSfc', xyz_DPTempDtFlux(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DExnerDtSfc', xyz_DExnerDtFlux(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelXDtSfc',  pyz_DVelXDtFlux (1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtSfc',  xqz_DVelYDtFlux (1:nx,1:ny,1:nz))
    do s = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(s))//'DtSfc', &
        & xyzf_DQMixDtFlux(1:nx,1:ny,1:nz,s))
    end do

    call HistoryAutoPut(TimeN, 'PTempSfcFlux', xyr_PTempFlux(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerSfcFlux', xyr_ExnerFlux(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'VelXSfcFlux',  pyr_VelXFlux (1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'VelYSfcFlux',  xqr_VelYFlux (1:nx,1:ny,1:nz))
    do s = 1, ncmax
      call HistoryAutoPut(TimeN, trim(SpcWetSymbol(s))//'SfcFlux', &
        & xyrf_QMixFlux(1:nx,1:ny,1:nz,s))
    end do

    call HistoryAutoPut(TimeN, 'SfcHeatFlux', xyr_SfcHeatFlux(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'SfcXMomFlux', pyr_SfcXMomFlux(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'SfcXMomFlux', xqr_SfcYMomFlux(1:nx,1:ny,1:nz))
    do s = 1, ncmax
      call HistoryAutoPut(TimeN, trim(SpcWetSymbol(s))//'SfcMassFlux', &
        & xyz_DensBZ(1:nx,1:ny,1:nz) * xyrf_QMixFlux(1:nx,1:ny,1:nz,s))
    end do

  end subroutine Surfaceflux_Const_forcing
  
end module Surfaceflux_Const
