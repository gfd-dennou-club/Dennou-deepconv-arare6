!= 簡単放射: とある高度を一様冷却
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: radiation_simple.f90,v 1.10 2014/03/04 04:49:41 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Radiation_Simple
  !
  ! 簡単放射: とある高度を一様冷却

  !モジュール読み込み
  !
  use dc_types, only: DP, STRING
  
  !暗黙の型宣言禁止
  !
  implicit none

  !private 属性
  !
  private

  !変数定義
  real(DP), save, allocatable, public :: xyz_DPTempDtRadVary(:,:,:)   
                                           !放射加熱項
  real(DP), save, allocatable, public :: xyz_DPTempDtRadConst(:,:,:)  
                                           !放射加熱項
  real(DP), save, allocatable, public :: xyz_ExnerRadVary(:,:,:)   
                                           !圧力方程式非断熱加熱項
  real(DP), save, allocatable, public :: xyz_ExnerRadConst(:,:,:) 
                                           !圧力方程式非断熱加熱項 

  real(DP), save               :: FactorDExnerDtRad = 1.0d0  
  character(STRING), parameter :: module_name = 'radiation_simple'
                                           ! モジュールの名称.
                                           ! Module name

  public Radiation_Simple_init
  public Radiation_HeatConst_forcing
  public Radiation_HeatVary_forcing

contains

!!!------------------------------------------------------------------------!!!
  subroutine Radiation_Simple_init
    !
    !NAMELIST から冷却率, 高度領域を設定.
    !

    ! モジュールの読み込み
    !
    use dc_types,          only : DP
    use dc_iounit,         only : FileOpen
    use dc_message,        only : MessageNotify
    use gtool_historyauto, only : HistoryAutoAddVariable
    use gridset,           only : imin,     & !x 方向の配列の下限
      &                           imax,     & !x 方向の配列の上限
      &                           jmin,     & !y 方向の配列の下限
      &                           jmax,     & !y 方向の配列の上限
      &                           kmin,     & !z 方向の配列の下限
      &                           kmax        !z 方向の配列の上限
    use axesset,           only : z_Z         !Z 座標軸(スカラー格子点)
    use constants,         only : DayTime     ! 1 日の長さ [s]
    use basicset,          only : xyz_ExnerBZ !エクスナー関数の基本場
    use DExnerDt,          only : xyz_DExnerDt_xyz
    use namelist_util,     only : namelist_filename
    
    !暗黙の型宣言禁止
    implicit none

    !入力変数
    real(DP) :: HeightUp   = 0.0d0  !放射強制を与える鉛直領域の上限
    real(DP) :: HeightDown = 0.0d0  !放射強制を与える鉛直領域の下限
    real(DP) :: RadHeatRate = 0.0d0 !一様放射加熱率 [K/day]
    integer  :: k                   !ループ変数
    integer  :: unit

    ! NAMELIST から情報を取得
    NAMELIST /radiation_simple_nml/ &
      & RadHeatRate, HeightUp, HeightDown, &
      & FactorDExnerDtRad

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=radiation_simple_nml)
    close(unit)

    allocate( xyz_DPTempDtRadVary(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_DPTempDtRadConst(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_ExnerRadVary(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_ExnerRadConst(imin:imax, jmin:jmax, kmin:kmax) )

    
    ! 温位の式の放射強制項. 
    ! 地表面から RadHeight で指定された高度までの間で一様放射冷却を与える. 
    !
    do k = kmin, kmax
      if ( z_Z(k) <= HeightDown  ) then
        xyz_DPTempDtRadConst(:,:,k) = 0.0d0 
      elseif( z_Z(k) >= HeightUp ) then
        xyz_DPTempDtRadConst(:,:,k) = 0.0d0 
      else
        xyz_DPTempDtRadConst(:,:,k) = RadHeatRate / DayTime / xyz_ExnerBZ(:,:,k)
      end if
    end do    
    xyz_ExnerRadConst = xyz_DExnerDt_xyz(xyz_DPTempDtRadConst) * FactorDExnerDtRad
    

    ! 温位の放射強制項.
    ! 地表面から HeightDown までは RadHeatRate で冷却. 
    ! HeightUp より上空は加熱率ゼロになるように加熱率を減少させる.
    ! Nakajima and Matsuno(1988),中島(1944)を参考にした
    !
    do k = kmin, kmax
      if ( z_Z(k) <= HeightDown ) then
        xyz_DPTempDtRadVary(:,:,k) = RadHeatRate / DayTime / xyz_ExnerBZ(:,:,k)
        
      elseif ( z_Z(k) > HeightDown .AND. z_Z(k) <= HeightUP) then
        xyz_DPTempDtRadVary(:,:,k) =  &
          & RadHeatRate / DayTime / xyz_ExnerBZ(:,:,k) &
          & * (HeightUP - z_Z(k)) / (HeightUP - HeightDown) 
        
      else if (z_Z(k) > HeightUP) then
        xyz_DPTempDtRadVary(:,:,k) = 0.0d0
      end if
    end do
    xyz_ExnerRadVary = xyz_DExnerDt_xyz(xyz_DPTempDtRadVary) * FactorDExnerDtRad

    ! Output
    !
    call MessageNotify( "M", &
      & module_name, "RadHeatRate = %f", d=(/RadHeatRate/))
    call MessageNotify( "M", &
      & module_name, "HeightUp = %f", d=(/HeightUP/))
    call MessageNotify( "M", &
      & module_name, "HeightDown= %f", d=(/HeightDown/))
    call MessageNotify( "M", &
      & module_name, "FactorDExnerDtRad= %f", d=(/ FactorDExnerDtRad /))

    ! ヒストリデータ定義
    ! 
    call HistoryAutoAddVariable(  &
      & varname='DPTempDtRad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Radiation term of potential temperature', &
      & units='K.s-1',    &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='ExnerRad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Radiation term of Exner function', &
      & units='K.s-1',    &
      & xtype='float')

  end subroutine Radiation_Simple_init


  subroutine Radiation_HeatConst_forcing(  &
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
    xyz_DPTempDt = xyz_DPTempDt + xyz_DPTempDtRadConst
    xyz_DExnerDt = xyz_DExnerDt + xyz_ExnerRadConst

    ! ファイル出力
    !    
    call HistoryAutoPut(TimeN, 'DPTempDtRad', xyz_DPTempDtRadConst(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerRad', xyz_ExnerRadConst(1:nx,1:ny,1:nz))   

  end subroutine Radiation_HeatConst_forcing


  subroutine Radiation_HeatVary_forcing(   &
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
    xyz_DPTempDt = xyz_DPTempDt + xyz_DPTempDtRadVary
    xyz_DExnerDt = xyz_DExnerDt + xyz_ExnerRadVary
    
    ! ファイル出力
    !    
    call HistoryAutoPut(TimeN, 'DPTempDtRad', xyz_DPTempDtRadVary(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerRad', xyz_ExnerRadVary(1:nx,1:ny,1:nz))
    
  end subroutine Radiation_HeatVary_forcing


end module Radiation_Simple
