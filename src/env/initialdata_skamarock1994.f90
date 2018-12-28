!= Module initialdata_Skamarock1994
!
! Authors::   SUGIYAMA Koichiro, ODAKA Masatsugu
! Version::   $Id: initialdata_skamarock1994.f90,v 1.1 2014/07/10 01:15:50 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module initialdata_Skamarock1994
  !
  ! Skamarock and Klemp (1994) を模した初期場を作成するためのモジュール

  !暗黙の型宣言禁止
  implicit none

  !初期化だけ公開
  public initialdata_SK1994_basic
  public initialdata_SK1994_disturbXZ
  public initialdata_SK1994_disturbYZ

contains

!!------------------------------------------------------------------------------!!!
  subroutine initialdata_SK1994_basic( z_Temp, z_Press )
    !
    !== 概要
    !  * Skamarock & Klemp (1994) で使われた初期値を再現する.
    !    * 浮力振動数一定 (N = 10^-2 1/sec)
    !    * N^2 = \frac{g}{\theta}{\DP{\theta}{z}}
    !

    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: kmin,          &!配列の Z 方向の下限
      &                   kmax,          &!配列の Z 方向の上限
      &                   nz              !配列の Z 方向の数
    use axesset, only:    z_Z, dz         !高度
    use constants, only:  GasRDry,       &!乾燥成分の定圧比熱
      &                   PressSfc,      &!地表面圧力
      &                   TempSfc,       &!地表面圧力
      &                   CpDry,         &!乾燥断熱減率
      &                   Grav            !重力
    
    !暗黙の型宣言禁止
    !
    implicit none

    !変数の定義
    !
    real(DP), intent(out) :: z_Press(kmin:kmax)       !圧力
    real(DP), intent(out) :: z_Temp(kmin:kmax)        !温度
    real(DP)              :: z_PTemp(kmin:kmax)
!    real(DP), parameter   :: PTemp0 = 300.0d0
!    real(DP), parameter   :: Press0 = 1.0d5
    real(DP), parameter   :: NN     = 1.0d-4
    integer               :: k

    ! 初期化
    !
    z_Press  = 0.0d0
    z_Temp   = 0.0d0
    z_PTemp  = 0.0d0

    ! 大気安定度をから温位を計算する
    ! 
    z_PTemp(1) = TempSfc
    
    do k = 2, nz
      z_PTemp(k) = TempSfc * dexp( NN * ( z_Z(k) - z_Z(1) ) / Grav )
    end do

    !----------------------------------------------
    ! 圧力と温度の計算
    !    
    z_Press(1) = PressSfc
    z_Temp(1)  = TempSfc
    
    do k = 2, nz
      z_Press(k) = z_Press(k-1) - ( Grav * z_Press(k-1) * dz ) / ( GasRDry * z_Temp(k-1) )
      z_Temp(k)  = z_PTemp(k) * ((z_Press(k) / PressSfc) ** (GasRDry / CpDry))
    end do
    
  end subroutine Initialdata_SK1994_basic


  subroutine initialdata_SK1994_disturbXZ(DelMax, Xc, Xr, xyz_Var)
    !
    ! 温位擾乱は θ(x,z) = A sin (πz/H) / [1 + (x-x_c)^2/a^2]
    ! A = 1.0e-2 K, x_c = 100 km, a = 5 km

    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: imin, imax,    &!配列の X 方向の上限・下限
      &                   jmin, jmax,    &!配列の Y 方向の上限・下限
      &                   kmin, kmax      !配列の Z 方向の上限・上限
    use axesset, only:    x_X, z_Z, ZMax     
    
    !暗黙の型宣言禁止
    !
    implicit none
    
    !変数の定義
    !
    real(DP), intent(in)   :: DelMax, Xc, Xr
    real(DP), intent(out)  :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), parameter    :: PI = 3.141592d0
    integer                :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_Var(i,j,k) = DelMax * sin(PI * z_Z(k) / Zmax ) / (1.0d0 + (x_X(i) - Xc)**2 / Xr**2 )
        end do
      end do
    end do
    
  end subroutine initialdata_SK1994_disturbXZ


  subroutine initialdata_SK1994_disturbYZ(DelMax, Yc, Yr, xyz_Var)
    !
    ! 温位擾乱は θ(x,z) = A sin (πz/H) / [1 + (x-x_c)^2/a^2]
    ! A = 1.0e-2 K, x_c = 100 km, a = 5 km

    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: imin, imax,    &!配列の X 方向の上限・下限
      &                   jmin, jmax,    &!配列の Y 方向の上限・下限
      &                   kmin, kmax      !配列の Z 方向の上限・上限
    use axesset, only:    y_Y, z_Z, ZMax     
    
    !暗黙の型宣言禁止
    !
    implicit none
    
    !変数の定義
    !
    real(DP), intent(in)   :: DelMax, Yc, Yr
    real(DP), intent(out)  :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), parameter    :: PI = 3.141592d0
    integer                :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_Var(i,j,k) = DelMax * sin(PI * z_Z(k) / Zmax ) / (1.0d0 + (y_Y(j) - Yc)**2 / Yr**2 )
        end do
      end do
    end do
    
  end subroutine initialdata_SK1994_disturbYZ

  
end module Initialdata_Skamarock1994
