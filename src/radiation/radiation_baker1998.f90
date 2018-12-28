!= Module Radiation
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: radiation_baker1998.f90,v 1.3 2014/03/04 05:55:05 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Radiation_baker1998
  !
  ! Baker et al (1998) を模した放射強制を与えるためのモジュール
  !

  ! モジュール読み込み
  !
  use dc_types, only: DP, STRING
  
  ! 暗黙の型宣言禁止
  !
  implicit none

  ! private 属性
  !
  private

  !変数定義
  !
  real(DP), save, allocatable, public :: xyz_DPTempDtRad(:,:,:)  !放射加熱が存在する領域
  real(DP), save, allocatable, public :: xyz_ExnerRad(:,:,:)  !放射加熱が存在する領域

  character(STRING), parameter :: module_name = 'radiation_baker1998'
                               ! モジュールの名称.
                               ! Module name

  public Radiation_baker1998_init
  public Radiation_baker1998_forcing

contains

!!!------------------------------------------------------------------------!!!
  subroutine Radiation_baker1998_init
    !
    ! 初期化ルーチン: BAKER1998 で用いられている熱強制を用いる
    !

    ! モジュール呼び出し
    !
    use dc_types,          only : DP
    use gtool_historyauto, only : HistoryAutoAddVariable
    use gridset,           only : imin,     & !x 方向の配列の下限
      &                           imax,     & !x 方向の配列の上限
      &                           jmin,     & !y 方向の配列の下限
      &                           jmax,     & !y 方向の配列の上限
      &                           kmin,     & !z 方向の配列の下限
      &                           kmax        !z 方向の配列の上限
    use axesset,           only : z_Z         !Z 座標軸(スカラー格子点)
    use constants,         only : CpDry       !比熱
    use basicset,          only : xyz_ExnerBZ, & !エクスナー関数の基本場
      &                           xyz_DensBZ     !基本場の密度
    use DExnerDt,          only : xyz_DExnerDt_xyz
    
    ! 暗黙の型宣言禁止
    !
    implicit none

    !内部変数
    !
    integer  :: k             !ループ変数

    allocate( xyz_DPTempDtRad(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_ExnerRad(imin:imax, jmin:jmax, kmin:kmax) )
    
    ! 温位の放射強制項.
    ! BAKER1998 で用いられている熱強制のプロファイルを用いる
    !
    do k = kmin, kmax
      xyz_DPTempDtRad(:,:,k) = cal_Qsub( z_Z(k) ) / CpDry / xyz_DensBZ(:,:,k) / xyz_ExnerBZ(:,:,k)
    end do
    xyz_ExnerRad = xyz_DExnerDt_xyz( xyz_DPTempDtRad ) 

    ! Output
    !
    call HistoryAutoAddVariable(  &
      & varname='DPTempDtRad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Radiation term of potential temperature', &
      & units='K.s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='ExnerRad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Radiation term of Exner function', &
      & units='K.s-1',    &
      & xtype='double')

  end subroutine Radiation_baker1998_init


  subroutine Radiation_baker1998_forcing(  &
    & xyz_DPTempDt, xyz_DExnerDt           & !(inout)
    & )

    ! モジュール読み込み
    !
    use dc_types,          only : DP
    use gtool_historyauto, only : HistoryAutoPut
    use timeset,           only : TimeN
    use gridset,           only : imin,     & !x 方向の配列の下限
      &                           imax,     & !x 方向の配列の上限
      &                           jmin,     & !y 方向の配列の下限
      &                           jmax,     & !y 方向の配列の上限
      &                           kmin,     & !z 方向の配列の下限
      &                           kmax,     & !z 方向の配列の上限
      &                           nx, ny, nz

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 入出力変数
    !
    real(DP), intent(inout) :: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDt(imin:imax,jmin:jmax,kmin:kmax)

    ! tendency の更新
    !
    xyz_DPTempDt = xyz_DPTempDt + xyz_DPTempDtRad
    xyz_DExnerDt = xyz_DExnerDt + xyz_ExnerRad
    
    ! ファイル出力
    !    
    call HistoryAutoPut(TimeN, 'DPTempDtRad', xyz_DPTempDtRad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerRad', xyz_ExnerRad(1:nx,1:ny,1:nz))
    
  end subroutine Radiation_baker1998_forcing
  

  function cal_Qsub( z0 )
    ! 
    ! BS1989 で用いられていた熱強制の式
    ! 同じルーチンが surfaceflux_baker1998.f90 にも含まれているので注意すること.

    ! モジュール読み込み
    !
    use dc_types,          only : DP

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 変数定義
    !
    real(DP), intent(in) :: z0
    real(DP)             :: cal_Qsub
    real(DP), parameter  :: z_Upper = 6.70d4      ! 上部頂点高度
    real(DP), parameter  :: z_Lower = 2.70d4      ! 下部頂点高度
    real(DP), parameter  :: c_Upper = 2.70d-2     ! 第 1 係数
    real(DP), parameter  :: c_Lower = 3.60d-3     ! 第 2 係数
    real(DP), parameter  :: s_Upper = 7.50d3      ! 第 1 標準偏差
    real(DP), parameter  :: s_Lower = 1.30d4      ! 第 2 標準偏差

    cal_Qsub =                                                                      &
      & (                                                                           &
      &     c_Lower * exp( - ( z0 - z_Lower )**2.0d0 / ( 2.0d0 * s_Lower **2.0d0) ) &
      &   + c_Upper * exp( - ( z0 - z_Upper )**2.0d0 / ( 2.0d0 * s_Upper **2.0d0) ) &
      &  ) 
    
    
  end function cal_Qsub

end module Radiation_baker1998
