!= 微分演算モジュール (2 次精度中央差分)
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: differentiate_center2.f90,v 1.4 2007-04-11 11:59:58 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module differentiate_center2
  !
  ! 2 次精度の微分演算を行うための関数群をまとめたパッケージ型モジュール
  ! 水平 Arakawa-C, 鉛直 Lorentz グリッドでの微分計算
  !

  !暗黙の型宣言禁止
  implicit none

  !private 属性
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

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function xa_dx_pa( pa_var ) 
    !
    ! x 微分: 半整数格子点上の微分
    !
    
    !モジュール読み込み
    use dc_types,  only : DP
    use gridset,   only : imin, imax, kmin, kmax
                                !配列サイズ
    use axesset,   only : dx    !格子点間隔
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in) :: pa_var(imin:imax, kmin:kmax)
                                     !微分演算の対象となる変数
    real(DP)             :: xa_dx_pa(imin:imax, kmin:kmax)
                                     !微分値
    integer             :: ix
    
    ! 1 階微分の計算
    !
    do ix = imin+1, imax
      xa_dx_pa(ix,:) = ( pa_Var(ix,:) - pa_Var(ix-1,:) ) / dx
    end do
    
    xa_dx_pa(imin,:) = xa_dx_pa(imin+1,:)
    
  end function xa_dx_pa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function pa_dx_xa( xa_var ) 
    !
    ! x 微分: 整数格子点上の微分
    !

    !モジュール読み込み
    use dc_types,  only : DP
    use gridset,   only : imin, imax, kmin, kmax
                                !配列サイズ
    use axesset,   only : dx    !格子点間隔

    !暗黙の型宣言禁止
    implicit none

    !変数定義    
    real(DP), intent(in) :: xa_var(imin:imax, kmin:kmax)
                                     !微分演算の対象となる変数
    real(DP)             :: pa_dx_xa(imin:imax, kmin:kmax)
                                     !微分値    
    integer             :: ix

    ! 1 階微分の計算
    !
    do ix = imin, imax-1
      pa_dx_xa(ix,:) = (xa_Var(ix+1,:) - xa_Var(ix,:))/dx
    end do
    
    pa_dx_xa(imax,:) = pa_dx_xa(imax-1,:)

  end function pa_dx_xa


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function az_dz_ar(ar_var) 
    !
    ! z 微分: 半整数格子点上の微分
    !

    !モジュール読み込み
    use dc_types,  only : DP
    use gridset,   only : imin, imax, kmin, kmax
                                !配列サイズ
    use axesset,   only : dz    !格子点間隔

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in) :: ar_var(imin:imax, kmin:kmax)
                                     !微分演算の対象となる変数
    real(DP)             :: az_dz_ar(imin:imax, kmin:kmax)
                                     !微分値
    integer             :: kz
    
    ! 1 階微分の計算
    !
    do kz = kmin+1, kmax
      az_dz_ar(:,kz) = ( ar_Var(:,kz) - ar_Var(:,kz-1) ) / dz
    end do
    
    az_dz_ar(:,kmin) = az_dz_ar(:,kmin+1)
    
  end function az_dz_ar

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function ar_dz_az( az_var ) 
    !
    ! z 微分: 整数格子点上の微分
    !

    !モジュール読み込み
    use dc_types,  only : DP
    use gridset,   only : imin, imax, kmin, kmax
                                !配列サイズ
    use axesset,   only : dz    !格子点間隔

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in) :: az_var(imin:imax, kmin:kmax)
                                     !微分演算の対象となる変数
    real(DP)             :: ar_dz_az(imin:imax, kmin:kmax)
                                     !微分値
    integer             :: kz

    do kz = kmin, kmax-1
      ar_dz_az(:,kz) = ( az_Var(:,kz+1) - az_Var(:,kz) ) / dz
    end do
    
    ar_dz_az(:,kmax) = ar_dz_az(:,kmax-1)

  end function ar_dz_az

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
end module differentiate_center2
