!= 微分演算モジュール (4 次精度中央差分)
!
! Authors::   KITAMORI Taichi, SUGIYAMA Ko-ichiro
! Version::   $Id: differentiate_center4.f90,v 1.3 2006-09-21 03:01:02 odakker Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module differentiate_center4
  !
  !4 次精度の微分演算を行うための関数群をまとめたパッケージ型モジュール
  !水平 Arakawa-C, 鉛直 Lorentz グリッドでの微分計算
  !
  ! 差分式は以下の通り
  !
  !    du/dx(i) = 9/8*[u(i)-u(i-1)]/dx - 1/24*[u(i+1)-u(i-2)]/dx
  !

  !暗黙の型宣言禁止
  implicit none

  !属性の指定
  private 
  
  !公開要素
  public xz_dx_pz, xr_dx_pr
  public pz_dx_xz, pr_dx_xr
  public xz_dz_xr, pz_dz_pr
  public xr_dz_xz, pr_dz_pz

  !--------------------------------------------
  interface xz_dx_pz
    module procedure xa_dx_pa
  end interface xz_dx_pz

  interface xr_dx_pr
    module procedure xa_dx_pa
  end interface xr_dx_pr
  !--------------------------------------------
  interface pz_dx_xz
    module procedure pa_dx_xa
  end interface pz_dx_xz
  
  interface pr_dx_xr
    module procedure pa_dx_xa
  end interface pr_dx_xr
  !--------------------------------------------
  interface xz_dz_xr
    module procedure az_dz_ar
  end interface xz_dz_xr

  interface pz_dz_pr
    module procedure az_dz_ar
  end interface pz_dz_pr
  !--------------------------------------------
  interface xr_dz_xz
    module procedure ar_dz_az
  end interface xr_dz_xz

  interface pr_dz_pz
    module procedure ar_dz_az
  end interface pr_dz_pz
  !--------------------------------------------  

contains 

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function xa_dx_pa( pa_var ) 
    !
    ! x 方向: 半整数格子上の微分演算
    !

    !モジュール呼び出し
    use dc_types,        only : DP
    use gridset,         only : imin, imax, kmin, kmax
                                      !領域サイズ
    use axesset,         only : dx    !格子間隔
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in) :: pa_var(imin:imax, kmin:kmax)
                                     !微分演算の対象となる変数
    real(DP)             :: xa_dx_pa(imin:imax, kmin:kmax)
                                     !微分値
    integer             :: ix                          

    do ix = imin+2, imax-1
      xa_dx_pa(ix,:) = &
        &   ( ( pa_Var(ix,:)   - pa_Var(ix-1,:) ) * 9.0d0 / ( 8.0d0 * dx ) ) &
        & - ( ( pa_Var(ix+1,:) - pa_Var(ix-2,:) ) / ( 24.0d0 * dx ) )
    end do
    
    xa_dx_pa(imin:imin+1,:) = 1.0d10
    xa_dx_pa(imax,:) = 1.0d10

  end function xa_dx_pa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function pa_dx_xa( xa_var ) 
    !
    ! x 方向: 整数格子点上での微分演算
    !
   
    !モジュール呼び出し
    use dc_types,        only : DP
    use gridset,         only : imin, imax, kmin, kmax
                                      !領域サイズ
    use axesset,         only : dx    !格子間隔

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in) :: xa_var(imin:imax, kmin:kmax)
                                     !微分演算の対象となる変数
    real(DP)             :: pa_dx_xa(imin:imax, kmin:kmax)
                                     !微分値
    integer             :: ix                          
    
    do ix = imin+1, imax-2
      pa_dx_xa(ix,:) = &
        &   ( ( xa_Var(ix+1,:) - xa_Var(ix,  :) ) * 9.0d0 / (8.0d0* dx ) ) &
        & - ( ( xa_Var(ix+2,:) - xa_Var(ix-1,:) ) / ( 24.0d0 * dx ) )
    end do
    
    pa_dx_xa(imin,:) = 1.0d10
    pa_dx_xa(imax-1:imax,:) = 1.0d10

  end function pa_dx_xa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function az_dz_ar( ar_var ) 
    !
    ! z 方向: 半整数格子点上での微分演算
    !

    !モジュール呼び出し
    use dc_types,        only : DP
    use gridset,         only : imin, imax, kmin, kmax
                                      !領域サイズ
    use axesset,         only : dz    !格子間隔

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in) :: ar_var(imin:imax, kmin:kmax)
                                     !微分演算の対象となる変数
    real(DP)             :: az_dz_ar(imin:imax, kmin:kmax)
                                     !微分値
    integer             :: kz
    
    do kz = kmin+2, kmax-1
      az_dz_ar(:,kz) = &
        &   ( ( ar_Var(:,kz)   - ar_Var(:,kz-1) ) * 9.0d0 / ( 8.0d0 * dz ) ) &
        & - ( ( ar_Var(:,kz+1) - ar_Var(:,kz-2) ) / ( 24.0d0 * dz ) )
    end do
    
    az_dz_ar(:,kmin:kmin+1) = 1.0d10
    az_dz_ar(:,kmax) = 1.0d10
    
  end function az_dz_ar
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function ar_dz_az( az_var ) 
    !
    ! z 方向: 整数格子点上での微分演算
    !

    !モジュール呼び出し
    use dc_types,        only : DP
    use gridset,         only : imin, imax, kmin, kmax
                                      !領域サイズ
    use axesset,         only : dz    !格子間隔

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in) :: az_var(imin:imax, kmin:kmax)
                                     !微分演算の対象となる変数    
    real(DP)             :: ar_dz_az(imin:imax, kmin:kmax)
                                     !微分値
    integer             :: kz

    do kz = kmin+1, kmax-2
      ar_dz_az(:,kz) = &
        &   ( ( az_Var(:,kz+1) - az_Var(:,kz)   ) * 9.0d0 / ( 8.0d0 * dz ) ) &
        & - ( ( az_Var(:,kz+2) - az_Var(:,kz-1) ) / ( 24.0d0 * dz ) )
    end do
    
    ar_dz_az(:,kmin) = 1.0d10
    ar_dz_az(:,kmax-1:kmax) = 1.0d10
    
  end function ar_dz_az

 
end module differentiate_center4
