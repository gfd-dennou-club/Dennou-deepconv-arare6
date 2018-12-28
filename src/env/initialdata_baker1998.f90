!= Module initialdata_baker1998
!
! Authors::   SUGIYAMA Koichiro, ODAKA Masatsugu
! Version::   $Id: initialdata_baker1998.f90,v 1.2 2014/03/04 05:55:04 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module initialdata_baker1998
  !
  ! Baker et al. (1998) を模した初期場を作成するためのモジュール

  !暗黙の型宣言禁止
  implicit none

  !デフォルトは private
  private

  !初期化だけ公開
  public initialdata_baker1998_basic

contains

!!------------------------------------------------------------------------------!!!
  subroutine initialdata_baker1998_basic( z_Temp, z_Press )
    !
    !== 概要
    !  * Baker and Shubert (1989) で使われた初期値を再現する.
    !    * 大気安定度から温度・圧力の基本場を作る
    !      * 温位のデータは安定度で示されているので, それを温位で焼き直す
    !      * 温位から, 基本場の温度・圧力を計算する. 
    !
    
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

    ! 初期化
    !
    z_Press  = 0.0d0
    z_Temp   = 0.0d0

    ! 大気安定度を与える.
    !   BS1998 の表 1 を目で読み取り, それを直線近似している. 
    ! 
    do k = kmin, kmax
      if ( r_Z(k) < 45.0d3 ) then 
        r_DPTempDz(k) = r_DPTempDz_45km + 0.30d-6 * ( r_Z(k) - 45.0d3 )
      elseif ( 45.0d3 <= r_Z(k) .AND. r_Z(k) < 48.0d3 ) then 
        r_DPTempDz(k) = r_DPTempDz_45km - 1.30d-6 * ( r_Z(k) - 45.0d3 )
      elseif ( 48.0d3 <= r_Z(k) .AND. r_Z(k) < 55.0d3 ) then 
        r_DPTempDz(k) = r_DPTempDz_48km
      elseif ( 55.0d3 <= r_Z(k) .AND. r_Z(k) < 60.0d3 ) then 
        r_DPTempDz(k) = r_DPTempDz_55km + 1.67d-6 * ( r_Z(k) - 55.0d3 )
      elseif ( 60.0d3 <= r_Z(k) ) then 
        r_DPTempDz(k) = r_DPTempDz_60km - 1.20d-6 * ( r_Z(k) - 60.0d3 )
      end if
    end do

    ! 温位勾配 DPTempDz を積分し温度場 z_PTemp を計算
    ! 傾きは半整数格子点の値を用いる. 
    !
    z_PTemp(nz) = TempTop - r_DPTempDz(1) * dz * 0.5d0
    do k = nz-1, 1, -1
      z_PTemp(k) = z_PTemp(k+1) - r_DPTempDz(k) * dz
    end do

    !----------------------------------------------
    ! 圧力と温度の計算
    !    
    z_Press(nz) = 2.20d4
    z_Temp(nz)  = 268.0d0
!    z_Press(nz) = PressTop + (Grav * PressTop * dz * 5.0d-1) / (GasRDry * TempTop)
!    z_Temp(nz)  = z_PTemp(nz) * (z_Press(nz) / PressBasis) ** (GasRDry / CpDry)

    do k = nz-1, 1, -1
      z_Press(k) = z_Press(k+1) + (Grav * z_Press(k+1) * dz) / (GasRDry * z_Temp(k+1))
      z_Temp(k)  = z_PTemp(k) * ((z_Press(k) / PressBasis) ** (GasRDry / CpDry))
    end do
       
  end subroutine Initialdata_baker1998_basic
  
end module Initialdata_baker1998
