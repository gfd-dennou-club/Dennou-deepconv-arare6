!= Module initialdata_smoothing
!
! Authors::   SUGIYAMA Koichiro, ODAKA Masatsugu
! Version::   $Id: initialdata_baker1998.f90,v 1.2 2014/03/04 05:55:04 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module initialdata_smoothing
  !
  ! 温位の鉛直プロファイルをスムージング.
  ! * サウンディングデータが「がたがた」している場合に対応するため.
  ! * 不安定成層な領域を中立成層に.

  !暗黙の型宣言禁止
  implicit none

  !デフォルトは private
  private

  !初期化だけ公開
  public initialdata_smoothing_basic

contains

!!------------------------------------------------------------------------------!!!
  subroutine initialdata_smoothing_basic( z_Temp, z_Press )
    !
    !== 概要
    ! 温位の鉛直プロファイルをスムージング.
    ! * サウンディングデータが「がたがた」している場合に対応するため.
    ! * 不安定成層な領域を中立成層に.
    
    !モジュール読み込み
    use dc_types,  only: STRING, DP
    use gridset,   only: kmin,          &!配列の Z 方向の下限
      &                  kmax,          &!配列の Z 方向の上限
      &                  nz              !配列の Z 方向の数
    use axesset,   only: r_Z,           &!高度
      &                  dz              !鉛直格子点間隔
    use constants, only: GasRDry,       &!乾燥成分の定圧比熱
      &                  TempTop,       &!地表面温度
!      &                  PressTop,      &!地表面圧力
      &                  PressBasis,    &!地表面圧力
      &                  CpDry,         &!乾燥断熱減率
      &                  Grav            !重力

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(out) :: z_Press(kmin:kmax)!圧力
    real(DP), intent(out) :: z_Temp(kmin:kmax) !温度
    real(DP)              :: r_DPTempDz(kmin:kmax)
    real(DP)              :: z_PTemp(kmin:kmax)
    real(DP), parameter   :: r_dPTempDz_45km = 3.9d-3
    real(DP), parameter   :: r_dPTempDz_48km = 0.0d0
    real(DP), parameter   :: r_dPTempDz_55km = 0.0d0
    real(DP), parameter   :: r_dPTempDz_60km = 8.35d-3
    integer               :: k


    ! 温位の傾きが負の場合には温位一定に.
    !
    do k = kmin+1, kmax
      if ( PTemp(k) <= PTemp(k-1) ) then 
        z_PTemp(k) = z_PTemp(k-1)
      end if
    end do

    !----------------------------------------------
    ! 圧力と温度の計算
    !    
    do k = kmin, kmax - 1
      z_Press(k) = z_Press(k+1) - (Grav * z_Press(k+1) * dz) / (GasRDry * z_Temp(k+1))
      z_Temp(k)  = z_PTemp(k) * ((z_Press(k) / PressBasis) ** (GasRDry / CpDry))
    end do
       
  end subroutine Initialdata_baker1998_basic
  
end module Initialdata_baker1998
