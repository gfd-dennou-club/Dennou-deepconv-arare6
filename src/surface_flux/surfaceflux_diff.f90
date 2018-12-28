!= Module HeatFlux
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: surfaceflux_diff.f90,v 1.16 2014/11/07 06:46:44 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Surfaceflux_diff
  !
  !下部境界でのフラックスの計算モジュール. 下部境界の拡散係数を決めうちする. 
  !
  
  !モジュール読み込み
  use dc_types, only: DP, STRING

  !暗黙の型宣言禁止
  implicit none

  !属性の指定
  private

  !関数を public に設定
  public surfaceflux_diff_init
  public surfaceflux_diff_forcing

  !変数定義
  real(DP), save  :: Kappa = 800.0d0               ! 下部境界での乱流拡散係数
  real(DP), save  :: FactorDExnerDtSurf = 1.0d0    ! Flag for diabatice heating term in pressure equation 
  character(STRING), parameter:: module_name = 'surfaceflux_diff'

contains
!!!------------------------------------------------------------------------!!!
  subroutine Surfaceflux_Diff_init
    !
    !NAMELIST から必要な情報を読み取り, 時間関連の変数の設定を行う. 
    !

    !モジュール読み込み
    !
    use dc_iounit,         only : FileOpen
    use dc_message,        only : MessageNotify
    use gtool_historyauto, only : HistoryAutoAddVariable
    use namelist_util,     only : namelist_filename
    use composition,       only : SpcWetSymbol
    use gridset,           only : ncmax

    !暗黙の型宣言禁止
    implicit none

    !内部変数
    integer    :: l, unit

    !---------------------------------------------------------------    
    ! NAMELIST から情報を取得
    !
    NAMELIST /surfaceflux_diff_nml/ &
      & Kappa,                      &
      & FactorDExnerDtSurf 

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=surfaceflux_diff_nml)
    close(unit)  

    call MessageNotify( "M", module_name, "Kappa = %f", d=(/Kappa/))
    call MessageNotify( 'M', module_name, "FactorDExnerDtSurf = %f", &
      &                  d=(/FactorDExnerDtSurf/) )

    call HistoryAutoAddVariable(  &
      & varname='DPTempDtFlux',         &
      & dims=(/'x','y','z','t'/), &
      & longname='surface flux of potential temperature', &
      & units='kg.kg-1.s-1',            &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='ExnerFlux',      &
      & dims=(/'x','y','z','t'/), &
      & longname='surface flux of Exner function', &
      & units='kg.kg-1.s-1',            &
      & xtype='double')

    do l = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtFlux', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='Surface Flux term of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio',  &
        & units='kg.kg-1.s-1',    &
        & xtype='double')
    end do

  end subroutine Surfaceflux_Diff_init


!!!------------------------------------------------------------------------!!!
  subroutine Surfaceflux_Diff_forcing( &
    &   xyz_PTemp, xyzf_QMix,          &
    &   xyz_DPTempDt, xyz_DExnerDt,    &
    &   xyzf_DQMixDt                   &
    & )
    ! 
    ! 下部境界からのフラックスによる温度の変化率を,
    ! バルク方法に基づいて計算する.
    !

    !モジュール読み込み
    !
    use dc_types,          only : DP
    use gtool_historyauto, only : HistoryAutoPut    
    use gridset,           only : imin,         & !x 方向の配列の下限
      &                           imax,         & !x 方向の配列の上限
      &                           jmin,         & !y 方向の配列の下限
      &                           jmax,         & !y 方向の配列の上限
      &                           kmin,         & !z 方向の配列の下限
      &                           kmax,         & !z 方向の配列の上限
      &                           nx, ny, nz, ncmax
    use axesset,           only : z_dz            !z 方向の格子点間隔
    use timeset,           only : TimeN
    use composition,       only : SpcWetSymbol, GasNum
    use DExnerDt,          only : xyz_DExnerDt_xyz_xyzf

    !暗黙の型宣言禁止
    !
    implicit none
    
    !変数定義
    !
    real(DP), intent(in)   :: xyz_PTemp(imin:imax,jmin:jmax,kmin:kmax)
                                           !温位の擾乱成分    
    real(DP), intent(in)   :: xyzf_QMix(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                           !温位の擾乱成分    
    real(DP), intent(inout):: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout):: xyz_DExnerDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout):: xyzf_DQMixDt(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP)               :: xyz_DPTempDt0(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)               :: xyz_DExnerDt0(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)               :: xyzf_DQMixDt0(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP)               :: xyz_Heatflux(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)               :: xyz_Exnerflux(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)               :: xyzf_QMixflux(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                            !地表面熱フラックス
    integer, parameter     :: kz = 1        !配列添字
    integer                :: l             !ループ変数

    ! 初期化
    !
    xyz_HeatFlux  = 0.0d0
    xyz_ExnerFlux = 0.0d0
    xyzf_QMixFlux = 0.0d0

    xyz_DPTempDt0 = xyz_DPTempDt
    xyz_DExnerDt0 = xyz_DExnerDt
    xyzf_DQMixDt0 = xyzf_DQMixDt

    ! 地表面熱フラックスによる加熱率を計算
    !   * 単位は K/s
    !   * エクスナー関数は基本場の値で代表させる.     
    !   * 格子点 xz では, 物理領域の最下端の添え字は kz = 1

    xyz_HeatFlux(:,:,kz) =                                  &
!     &  - Kappa * xyz_PTemp(:,:,kz) * xyz_ExnerBZ(:,:,kz)  &  !check
      &  - Kappa * xyz_PTemp(:,:,kz)                        &
      &    / ( ( z_dz(kz) * 5.0d-1 ) ** 2.0d0 )  

    do l = 1, GasNum
      xyzf_QMixFlux(:,:,kz,l) =                        &
        &     - Kappa * xyzf_QMix(:,:,kz,l)            &
        &        / ( ( z_dz(kz) * 5.0d-1 ) ** 2.0d0 )  
    end do
    
    xyz_DPTempDt = xyz_DPTempDt0 + xyz_Heatflux
    xyzf_DQMixDt = xyzf_DQMixDt0 + xyzf_QMixflux
    
    xyz_ExnerFlux = xyz_DExnerDt_xyz_xyzf( xyz_HeatFlux, xyzf_QMixFlux ) * FactorDExnerDtSurf
    xyz_DExnerDt = xyz_DExnerDt0 + xyz_ExnerFlux
    
    call HistoryAutoPut(TimeN, 'DPTempDtFlux', xyz_HeatFlux(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerFlux', xyz_ExnerFlux(1:nx,1:ny,1:nz))
    do l = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(l))//'DtFlux', xyzf_Qmixflux(1:nx,1:ny,1:nz,l))
    end do    

  end subroutine Surfaceflux_Diff_forcing
  
end module Surfaceflux_diff
